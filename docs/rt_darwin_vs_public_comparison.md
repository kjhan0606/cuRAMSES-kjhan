# DARWIN RT vs Public RT 비교 분석

## 요약

DARWIN private RT (`/gpfs/kjhan/Darwin/Cudaization/RAMSES-DARWIN-v2/`)는 public RT 대비 4가지 핵심 강화:
1. **GPU 가속**: CUDA 커널 + Fortran ISO_C_BINDING 인터페이스 (cooling + flux)
2. **고급 물리**: H2 자기차폐, subgrid HII 가열, 광전 가열, 우주선 이온화
3. **OpenMP 병렬화**: `patch/darwin_omp/`에 전체 RT 모듈 OMP 지원
4. **Morton key 최적화**: `patch/kjhan/`에서 O(1) 이웃 탐색

## 파일별 비교

### 동일한 파일 (1개)
- `rt_units.f90` — 변경 없음

### 변경된 파일 (18개)
| 파일 | 변경 라인 | 분류 |
|------|----------|------|
| rt_cooling_module.f90 | ~1361 | **대규모** — HII 가열, 우주선, H2 차폐 |
| rt_spectra.f90 | ~937 | **대규모** — 확장 스펙트럼 |
| coolrates_module.f90 | ~535 | cooling rate 테이블 구조 변경 |
| rt_output_hydro.f90 | ~426 | 출력 형식 |
| rt_init.f90 | ~288 | 이온 초기화, 가열 변수 |
| rt_init_xion.f90 | ~185 | 이온화 분율 초기화 |
| rt_flux_module.f90 | ~185 | flux 계산 업데이트 |
| rt_init_flow_fine.f90 | ~109 | IC 읽기 |
| rt_parameters.f90 | ~100 | 30+ 새 파라미터 |
| rt_interpol_hydro.f90 | ~58 | 보간 |
| rt_hydro_boundary.f90 | ~54 | 경계 조건 |
| rt_init_hydro.f90 | ~42 | 초기화 |
| rt_hydro_flag.f90 | ~22 | refinement 플래그 |
| rt_godunov_utils.f90 | ~18 | 유틸리티 |
| rt_boundana.f90 | ~12 | 경계 분석 |
| rt_condinit.f90 | ~4 | 포맷 |
| rt_hydro_commons.f90 | ~4 | 파라미터 정의 |

### DARWIN 전용 파일 (2개 — 핵심)

#### rt_cuda_interface.f90 (576 lines)
- ISO_C_BINDING 기반 Fortran-CUDA 인터페이스
- 주요 함수:
  - `rt_cuda_init_f()` — GPU 초기화 (cooling 테이블 + RT 파라미터 전송)
  - `rt_cuda_solve_cooling_f()` / `_async_f()` — cooling 커널 (동기/비동기)
  - `rt_cuda_compute_flux_f()` / `_async_f()` — flux 계산 (동기/비동기)
  - `rt_cuda_acquire_stream_f()` — 스트림 풀 관리
- 21개 cooling rate 테이블을 연속 배열로 flatten → GPU 디바이스 메모리 전송
- HLL eigenvalue 테이블 (101×101) 전송
- `MPI_Comm_split_type`로 multi-GPU affinity 지원

#### rt_cuda_kernels.cu (CUDA 소스)
- `rt_cuda_solve_cooling(...)` — 셀별 cooling ODE solver
- `rt_cuda_compute_flux(...)` — HLL Riemann solver
- BLOCK_SIZE=128, 스트림 풀 파이프라이닝
- `cuda_stream_pool.h` 의존 (이미 CUBE_HR5에 존재)

## DARWIN 추가 물리

### H2 자기차폐 (Self-Shielding)
```fortran
real(dp),dimension(1:NGROUPS)::ssh2 = 1d0   ! 차폐 인자
real(dp),dimension(1:NGROUPS)::isLW = 0d0   ! Lyman-Werner 밴드 소속
```
- 칼럼 밀도 기반 H2 해리 차폐
- Nickerson+ 2018 모델

### Subgrid HII 가열
```fortran
integer::heat_unresolved_HII=0  ! 0=off, 1/2=subgrid 모델
```
- `heat_unresolved_HII_regions(ilevel)` — 미분해 이온화 영역 가열
- NENER 슬롯 또는 passive scalar 사용

### 광전 가열 (Photoelectric Heating)
```fortran
real(dp),dimension(nGroups)::kappaAbs=0  ! 먼지 흡수 불투명도 [cm²/g]
real(dp),dimension(nGroups)::kappaSc=0   ! 먼지 산란 불투명도
integer::iPEH_group=0                    ! 광전 가열 광자 그룹
```

### 우주선 이온화 (Cosmic Ray)
```fortran
logical::cosmic_rays=.false.
real(dp),parameter::cosray_H2 = 7.525d-16  ! [s⁻¹] (Indriolo 2012)
real(dp),parameter::cosray_HI = 4.45d-16   ! [s⁻¹] (Indriolo 2015)
```

### AGN 복사
```fortran
logical::rt_AGN=.false.
```

## OpenMP 지원 (patch/darwin_omp/)

12개 RT 파일 전체에 OMP 지시자 추가:
- `!$omp parallel do private(...) schedule(dynamic)` 패턴
- `save` 배열 → thread-local 배열 전환 (`#ifdef _OPENMP`)
- `rt_godunov_fine.f90`, `rt_cooling_module.f90` 등

## Morton Key 최적화 (patch/kjhan/)

3개 파일:
- `rt_godunov_fine.f90` — `son(nbor())` → `hash_nbor()`, `hash_child_from_cell()`
- `rt_hydro_boundary.f90` — 경계 조건 이웃 탐색
- `rt_hydro_flag.f90` — refinement 플래그

## 포팅 우선순위

| 구성요소 | 라인 | 소요 | 효과 | 위험 | 추천 |
|---------|------|------|------|------|------|
| CUDA 가속 (interface+kernels) | 850 | 2-3일 | 매우 높음 (10-20×) | 낮음 | **YES** |
| H2 자기차폐 | 150 | 1일 | 높음 | 낮음 | **YES** |
| Subgrid HII 가열 | 200 | 2일 | 중간 | 중간 | **YES** |
| 광전 가열 | 150 | 1.5일 | 높음 | 낮음 | **YES** |
| 우주선 이온화 | 100 | 1일 | 중간 | 낮음 | 선택 |
| AGN 복사 | 100 | 1.5일 | 낮음 | 중간 | 나중에 |
| Morton key RT | 50 | 0.5일 | 중간 | 낮음 | 자동 (이미 있음) |
| OpenMP RT | 100 | 1일 | 중간 | 낮음 | 자동 |

**추천 총 소요: 7-8일** (Priority A+B: GPU + 핵심 물리)

## 소스 위치

- DARWIN GPU: `/gpfs/kjhan/Darwin/Cudaization/RAMSES-DARWIN-v2/rt/rt_cuda_interface.f90`
- DARWIN 물리: `/gpfs/kjhan/Darwin/Cudaization/RAMSES-DARWIN-v2/patch/darwin/rt_*.f90`
- DARWIN OMP: `/gpfs/kjhan/Darwin/Cudaization/RAMSES-DARWIN-v2/patch/darwin_omp/rt_*.f90`
- DARWIN kjhan: `/gpfs/kjhan/Darwin/Cudaization/RAMSES-DARWIN-v2/patch/kjhan/rt_*.f90`
- Public RT: `/gpfs/kjhan/Hydro/CUBE_HR5/code_cube/rt/`
