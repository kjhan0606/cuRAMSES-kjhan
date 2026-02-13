# Zoom-In Simulation Setup Guide (2026-02-14)

## 1. 개요

1 노드 (64 코어, 1 TB 메모리)에서 별 탄생이 가능한 zoom-in 우주론 시뮬레이션을
안정적으로 실행하기 위한 설정 과정을 기록한다.

### 목표
- Physical resolution: ~1 kpc (별 탄생 가시화)
- 1 노드에서 안정 실행 (메모리 < 1 TB)
- Zoom 지역에만 AMR 활성화, 배경 지역 AMR 비활성화

### 최종 결과
- Main step 15+까지 안정 실행
- 별 29개 생성 (Level 12-14)
- 메모리: 35.9 GB/proc (12 프로세스, 총 ~430 GB)
- AMR Level 14 도달 (dx_phys ~ 0.9 kpc at z=0)

---

## 2. 시뮬레이션 설계

### 박스 및 해상도
| 항목 | 값 |
|------|-----|
| Box size | 10 h^-1 Mpc |
| levelmin | 7 (base grid: 128^3) |
| levelmax | 14 (dx ~ 0.9 kpc physical at z=0) |
| IC levels | 7-11 (MUSIC로 생성) |
| Zoom region | 1.25 h^-1 Mpc cube (box 중심) |
| DM mass (finest) | ~8,759 h^-1 M_sun |

### 해상도 테이블
| Level | dx_comoving (h^-1 kpc) | dx_phys (kpc, z=0) | DM mass (h^-1 M_sun) |
|-------|----------------------|-------------------|---------------------|
| 7 | 78.1 | 115.4 | 배경 (heavy) |
| 8 | 39.1 | 57.7 | 배경 (heavy) |
| 9 | 19.5 | 28.9 | 배경 (heavy) |
| 10 | 9.77 | 14.4 | 배경 (heavy) |
| 11 | 4.88 | 7.22 | 8,759 (finest IC) |
| 12 | 2.44 | 3.61 | AMR only |
| 13 | 1.22 | 1.80 | AMR only |
| 14 | 0.61 | 0.90 | AMR only |

---

## 3. IC 생성 (MUSIC)

### MUSIC 설정 파일
- 경로: `test_ksection/zoomin_10Mpc_lv11.conf`
- MUSIC 바이너리: `/gpfs/kjhan/Darwin/Cudaization/MUSIC/build/MUSIC`

```
[setup]
boxlength       = 10          # Box size in Mpc/h
zstart          = 50          # Starting redshift
levelmin        = 7           # Base grid: 128^3
levelmax        = 11          # Zoom finest: 2048 effective
ref_center      = 0.5, 0.5, 0.5
ref_extent      = 0.1, 0.1, 0.1  # ~1.0 Mpc/h zoom region
baryons         = yes
use_2LPT        = yes

[cosmology]
Omega_m = 0.3111, Omega_L = 0.6889, Omega_b = 0.04
H0 = 67.66, sigma_8 = 0.8102, nspec = 0.9665

[random]
seed[7]  = 95064
seed[8]  = 31415
seed[9]  = 27182
seed[10] = 16180
seed[11] = 14142

[output]
format   = grafic2
filename = IC_zoomin_lv11
```

### 생성된 IC 구조
```
IC_zoomin_lv11/
  level_007/  (128x128x128, 전체 박스)
  level_008/  (66x66x66, zoom sub-region)
  level_009/  (96x96x96)
  level_010/  (152x152x152)
  level_011/  (256x256x256, finest zoom)
```

각 level에 포함된 파일:
`ic_deltab, ic_velbx/y/z, ic_poscx/y/z, ic_velcx/y/z, ic_pvar_00001, ic_refmap`

---

## 4. 핵심 문제: 배경 AMR 폭발

### 문제 현상
`ivar_refine=0` (quasi-Lagrangian) 설정 시, 배경 지역의 모든 셀이
level 7 → 8 → 9 → ... 로 기하급수적으로 refine되어 메모리 초과 크래시 발생.

### 원인 분석
`poisson_refine()` (amr/flag_utils.f90)의 `.not. init` 분기:
```fortran
if(.not. init) then
   do i=1,ncell
      ok(i)=ok(i).or.(cpu_map2(ind_cell(i))==1)
   end do
```
`cpu_map2`는 `rho_fine`에서 `phi >= m_refine`로 설정되는데,
배경 셀에도 입자가 있어 `cpu_map2=1`이 설정될 수 있음.

### 시도한 해결 방법
| 시도 | 결과 |
|------|------|
| `sink_refine=.false.` | 불충분 — cpu_map2 문제 미해결 |
| `ngridtot=200M` (80M→200M) | 불충분 — grid 폭발 속도만 지연 |
| `ivar_refine=-1` (density 기반) | init 단계에서 크래시 — IC 구조와 불일치 |
| **`ivar_refine=11` + `ic_pvar_00006`** | **성공** |

---

## 5. 해결: Zoom Geometry Passive Scalar (HR5 방식)

### 원리
HR5 시뮬레이션에서 사용된 방법:
- 6번째 passive scalar (`ic_pvar_00006`)에 zoom geometry 정보 저장
- `ivar_refine=11` (NVAR=11: hydro 5 + passive 6)
- `var_cut_refine=0.01`: scalar > 0.01인 셀만 refine 허용

### Refinement 판단 로직
```
poisson_refine에서:
  init 단계:    uold(cell, 11) / uold(cell, 1) > 0.01  → refine
  runtime 단계: cpu_map2(cell) == 1                      → refine
                (cpu_map2는 rho_fine에서 mass_cut_refine 적용 후 설정)
```

- **init 단계**: `ic_pvar_00006`에서 읽은 scalar 값으로 줌 셀만 식별
- **runtime 단계**: `mass_cut_refine`이 heavy DM 입자를 필터링하여
  cpu_map2가 줌 지역에서만 1로 설정됨

### ic_pvar_00006 생성
스크립트: `test_ksection/create_pvar006.py`

```python
# Level 7 (base): ic_refmap에서 zoom mask 추출 (0.2% = 4096/2097152 셀)
# Levels 8-11 (zoom): 모든 셀 = 1.0 (전체가 zoom 영역)
```

| Level | 파일 크기 | Zoom 셀 비율 |
|-------|----------|-------------|
| 7 | 128^3 | 0.2% (4096 셀) |
| 8 | 66^3 | 100% |
| 9 | 96^3 | 100% |
| 10 | 152^3 | 100% |
| 11 | 256^3 | 100% |

### init_refine.f90의 핵심 로직
```fortran
! amr/init_refine.f90 line 31
if(ivar_refine==0) call init_refmap   ! ivar_refine=11이면 호출 안 됨
```
`ivar_refine=11`일 때 `init_refmap` 미호출 → `cpu_map2`는 0으로 유지 →
init 단계에서는 scalar 기준으로만 refinement 결정.

---

## 6. 최종 Namelist

파일: `test_ksection/run_zoomin_lv11/cosmo_zoomin_physics.nml`

### REFINE_PARAMS (핵심)
```fortran
&REFINE_PARAMS
m_refine=11*8.
ivar_refine=11              ! passive scalar 6 = zoom geometry
var_cut_refine=0.01         ! scalar > 0.01 → refine
interpol_var=1
interpol_type=0
sink_refine=.true.
mass_cut_refine=2.32831e-10 ! level 11 DM mass threshold
/
```

### 전체 파라미터 요약
| 파라미터 | 값 | 설명 |
|---------|-----|------|
| levelmin | 7 | Base grid 128^3 |
| levelmax | 14 | 최대 AMR (dx ~ 0.9 kpc) |
| ngridtot | 200,000,000 | 전체 grid 할당 |
| nparttot | 50,000,000 | 전체 입자 할당 |
| ivar_refine | 11 | Zoom geometry scalar |
| var_cut_refine | 0.01 | Zoom scalar threshold |
| mass_cut_refine | 2.32831e-10 | Heavy DM 필터 (lmin=11) |
| ordering | ksection | K-section 도메인 분할 |
| memory_balance | .true. | 메모리 기반 부하 분산 |
| nremap | 5 | 부하 분산 주기 |
| nsubcycle | 1,1,1,1,2 | Fine level subcycling |

---

## 7. 실행 결과

### 실행 환경
- 12 MPI 프로세스 (64 코어 중 12 사용)
- 실행 디렉토리: `test_ksection/run_zoomin_lv11/`
- 바이너리: `bin/` 에서 빌드된 ramses3d

### Main Step 15 시점 Grid 구조
```
Level  7:   262,144 grids (base, 전체 박스)
Level  8:    59,348 grids (zoom only)
Level  9:   326,602 grids
Level 10: 2,057,638 grids
Level 11: 10,695,246 grids
Level 12:   112,139 grids (AMR)
Level 13:    14,211 grids (AMR)
Level 14:       729 grids (AMR, 최대 레벨)
```

### 물리 결과
- **별 탄생**: 29개 (Level 12: 12, Level 13: 8, Level 14: 9)
- **에너지 보존**: econs = -2.74E-03
- **Scale factor**: a = 0.0232 (z ~ 42)

### 성능
- 메모리: 35.9 GB/proc (min/max: 35.5-35.9 GB)
- Grid 사용률: 7.8% of ngridmax
- 성능: ~5 us/pt, ~39 s/coarse step
- 총 실행 시간: 583초 (Main step 15)

### 안정성 확인
- Grid 폭발 없음: Level 8이 59K에서 안정 (이전: 2.1M으로 폭발)
- 메모리 증가 없음: 35.3 GB → 35.9 GB (15 step 동안 +0.6 GB)
- Morton mismatch: 0

---

## 8. mass_cut_refine 참고값

HR5에서 사용된 mass_cut_refine 값 (levelmin별):

| IC finest level | mass_cut_refine |
|----------------|-----------------|
| 8 | 1.19209e-07 |
| 9 | 1.49012e-08 |
| 10 | 1.86265e-09 |
| 11 | 2.32831e-10 |
| 12 | 2.91038e-11 |
| 13 | 3.63798e-12 |

---

## 9. 파일 목록

| 파일 | 설명 |
|------|------|
| `test_ksection/zoomin_10Mpc_lv11.conf` | MUSIC IC 생성 설정 |
| `test_ksection/IC_zoomin_lv11/` | 생성된 IC (levels 7-11) |
| `test_ksection/create_pvar006.py` | ic_pvar_00006 생성 스크립트 |
| `test_ksection/run_zoomin_lv11/cosmo_zoomin_physics.nml` | RAMSES namelist |
| `test_ksection/run_zoomin_lv11/run_zoomin_lv11.log` | 실행 로그 |

---

## 10. 문제 해결 요약

1. **시스템 크래시 진단**: 이전 실행이 시스템 비정상 재부팅으로 중단됨 (프로그램 버그 아님)
2. **시뮬레이션 규모 설계**: 1 노드에서 실행 가능한 box=10 Mpc/h, IC level 7-11 구성
3. **IC 재생성**: MUSIC로 ref_extent=0.1 (1 Mpc/h zoom region) IC 생성
4. **Grid 폭발 문제**: `ivar_refine=0`에서 배경 전체가 refine → OOM
5. **HR5 방식 적용**: `ivar_refine=11` + `ic_pvar_00006` (zoom geometry scalar)로 해결
6. **안정 실행 확인**: Main step 15+, 별 탄생 29개, 메모리 35.9 GB/proc
