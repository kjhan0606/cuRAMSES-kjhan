# RAMSES-RT 통합 계획

## 현재 상태

RAMSES-RT (Rosdahl et al. 2013) 코드가 `rt/` 디렉토리에 이미 존재 (7,142 lines, 20 files).
현재 빌드(`patch/Horizon5-master-2`)에는 포함되지 않음. `-DRT` 컴파일 플래그 + `rt=.true.` namelist 필요.

## 기존 RT 코드 구조

### 핵심 모듈 (rt/ 디렉토리)
| 파일 | 라인 | 역할 |
|------|------|------|
| `rt_parameters.f90` | 141 | RT 파라미터 정의 (nGroups, 감소 광속, M1 closure 등) |
| `rt_hydro_commons.f90` | 6 | `rtuold`, `rtunew` 배열 선언 |
| `rt_spectra.f90` | 1,749 | 항성 SED (spectral energy distribution) |
| `rt_cooling_module.f90` | 1,047 | RT cooling/heating + ionization |
| `coolrates_module.f90` | 528 | cooling rate 테이블 |
| `rt_init.f90` | 558 | RT 초기화 |
| `rt_godunov_fine.f90` | 517 | M1 hyperbolic solver (Godunov method) |
| `rt_flux_module.f90` | 418 | RT flux 계산 (GLF/HLL) |
| `rt_interpol_hydro.f90` | 175 | RT restriction (upload_fine) |
| `rt_init_flow_fine.f90` | 419 | RT IC 설정 |
| `rt_output_hydro.f90` | 290 | RT 출력 |
| `rt_hydro_boundary.f90` | 181 | RT 경계 조건 |

### 메인 루프 통합 (amr_step.f90)
`#ifdef RT` 블록 6개:
1. Ghost zone exchange: `make_virtual_fine_dp(rtuold)` + `rt_make_boundary_hydro`
2. Star RT feedback: `update_star_RT_feedback`
3. RT state 초기화: `rt_set_unew`
4. **RT 메인 루프** (`rt_step`): subcycling + godunov solver + cooling
5. UV rate 업데이트
6. RT 통계 출력

### MPI 통신
- RT 전용 MPI 호출 없음 — AMR 인프라(`make_virtual_fine_dp`, `make_virtual_reverse_dp`) 재사용
- 변수: `rtuold(1:ncell, 1:nrtvar)`, nrtvar = nGroups × (1+ndim)
- 이미 k-section 최적화된 가상 경계 교환을 자동으로 사용

### 기존 테스트 문제 (patch/rt/)
- `davis/`: 가스 부양 실험 (Davis+2014)
- `davis_cDecr/`: Davis + cooling/dust
- `iliev6_profile/`: 1/r² 밀도 프로파일 RT 비교 테스트
- `mol/`: 분자/금속 화학 + UV 의존 cooling

## 통합 단계

### Phase 1: 컴파일 통합
1. `bin/Makefile` 수정:
   - `DEFINES`에 `-DRT -DNGROUPS=$(NGROUPS) -DNIONS=$(NIONS)` 추가 (조건부)
   - `VPATH`에 `../rt` 추가
   - `MODOBJ`에 RT 모듈 추가: `rt_parameters.o`, `rt_hydro_commons.o`, `coolrates_module.o`, `rt_spectra.o`, `rt_cooling_module.f90`, `rt_flux_module.o`
   - `RTOBJ` 추가: `rt_init_hydro.o`, `rt_init_xion.o`, `rt_init.o`, `rt_init_flow_fine.o`, `rt_output_hydro.o`, `rt_godunov_fine.o`, `rt_interpol_hydro.o`, `rt_godunov_utils.o`, `rt_condinit.o`, `rt_hydro_flag.o`, `rt_hydro_boundary.o`, `rt_boundana.o`, `rt_read_hydro_params.o`, `rt_units.o`
2. `amr_step.jaehyun.f90`에 `#ifdef RT` 블록 이식 (원본 amr_step.f90 참조)
3. 컴파일 확인: `make RT=1 NGROUPS=3 NIONS=3`

### Phase 2: K-section 최적화
- RT는 이미 `make_virtual_fine_dp`/`make_virtual_reverse_dp` 사용 → k-section 자동 적용
- Bulk exchange 적용: `rtuold` nrtvar 컬럼을 한 번에 교환
- 예상 효과: RT subcycle당 통신 횟수 nrtvar→1로 감소

### Phase 3: GPU 가속
- `rt_godunov_fine`는 hydro godunov와 동일 구조 (sweep-based)
- `hydro_cuda_kernels.cu`와 같은 패턴으로 RT 커널 작성 가능
- RT flux 계산 (GLF/HLL)이 GPU에 적합 (cell-local 연산)

### Phase 4: 물리 테스트
1. Iliev Test 6 (1/r² 프로파일) — 해석해와 비교
2. Davis 가스 부양 테스트
3. Cosmological reionization 테스트 (z_reion ~ 6-10)
4. Star+RT feedback이 은하 형성에 미치는 영향 확인

### Phase 5: 프로덕션 설정
- Photon group 수 결정 (보통 3-5: HI ionizing, HeI ionizing, HeII ionizing, UV, IR)
- 감소 광속(reduced speed of light) 설정 (`rt_c_fraction`)
- Star particle RT emission 연동 (`rt_star=.true.`)
- Cooling module과 RT cooling 통합 검증

## 주요 파라미터 (namelist)
```
&RT_PARAMS
rt_advect=.true.
rt_star=.true.
rt_esc_frac=1.0
rt_courant_factor=0.8
rt_c_fraction=0.01    ! 감소 광속 (1% of c)
rt_otsa=.true.         ! On-the-spot approximation
rt_flux_scheme='glf'   ! GLF or HLL
rt_nsubcycle=5         ! RT subcycle per hydro step
/
```

## 예상 난이도
- Phase 1 (컴파일): 중간 — Makefile + amr_step 수정
- Phase 2 (K-section): 낮음 — 이미 인프라 있음
- Phase 3 (GPU): 높음 — 새 CUDA 커널 필요
- Phase 4-5 (테스트/프로덕션): 높음 — 물리 검증 필요

## 참고 논문
- Rosdahl et al. 2013, MNRAS 436, 2188 (RAMSES-RT)
- Rosdahl & Teyssier 2015, MNRAS 449, 4380 (subcycling)
- Levermore 1984, JQSRT 31, 149 (M1 closure)
