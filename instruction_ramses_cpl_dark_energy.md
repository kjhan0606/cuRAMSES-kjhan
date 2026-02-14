# RAMSES에 CPL 암흑에너지 구현 지시문

## 프로젝트 개요

RAMSES (AMR cosmological hydrodynamics code)에 CPL 파라미터화
`w(a) = w0 + wa*(1-a)`를 따르는 동적 암흑에너지(Dark Energy) 모형을 구현한다.

DE 섭동 처리는 음속(sound speed) 값에 따라 두 가지로 나뉜다:
- **Case A** (`c_s² = 0`): DE가 pressureless fluid처럼 클러스터링 → 실공간 유체 방정식
- **Case B** (`c_s² ≠ 0`): 준정적 근사 → screened Poisson (Helmholtz) 방정식을 multigrid로 풀기

**핵심 원칙:**
- `w0 = -1, wa = 0`이면 ΛCDM과 bit-wise 동일한 결과를 재현해야 함
- 기존 코드 수정을 최소화하고, 전처리기 지시문(`#ifdef`)으로 선택 가능하게 할 것
- 기존 RAMSES의 코딩 스타일(Fortran 90)을 준수할 것

---

## 단계 1: 코드 구조 파악

RAMSES 코드의 전체 구조를 파악해줘. 특히 다음을 분석해:

1. **Hubble function H(a) 정의 위치**: 어떤 파일, 어떤 서브루틴에서 H(a)를 계산하는지
2. **배경 우주론 파라미터**: `Omega_m`, `Omega_Lambda`, `H0` 등이 선언/초기화되는 위치
3. **Poisson solver**: multigrid 소스 파일 위치, 주요 서브루틴명
   - relaxation sweep (Gauss-Seidel) 루틴
   - restriction / prolongation 연산자
   - V-cycle 또는 FMG cycle 구조
   - 소스항(rho) 할당 부분
4. **입자 시간 적분기**: kick-drift-kick 루틴, 스텝 크기 계산
5. **바리온 수력학 솔버**: Godunov solver, 냉각 함수, UV 배경 등
6. **namelist 파라미터**: 읽기/파싱 루틴 위치
7. **초기조건 관련 루틴**: 성장인자 D(a) 계산, IC 생성
8. **AMR 구조**: 레벨별 격자 변수 선언, refinement criteria

각 파일의 핵심 서브루틴과 변수명을 정리해줘.

---

## 단계 2: 배경 우주론 수정 (공통, Case A/B 모두 필요)

RAMSES에 CPL 암흑에너지 배경 우주론을 구현해줘.

### 상태방정식
```
w(a) = w0 + wa * (1 - a)
```

### 수정 사항

#### 2-1. namelist 파라미터 추가
- `w0_de`: DE 상태방정식 현재값 (기본값 = -1.0)
- `wa_de`: DE 상태방정식 시간변화율 (기본값 = 0.0)
- `cs2_de`: DE 유효 음속 제곱 (기본값 = 1.0)
- `de_perturb`: DE 섭동 활성화 여부 (기본값 = .false.)

기본값 설정 시 `w0=-1, wa=0`이 되어 ΛCDM과 동일해야 한다.

#### 2-2. Hubble function 교체
기존:
```
H²(a) = H0² [ Omega_m * a^(-3) + Omega_Lambda ]
```
수정:
```
H²(a) = H0² [ Omega_m * a^(-3) + Omega_de * f_de(a) ]

f_de(a) = a^( -3*(1 + w0 + wa) ) * exp( -3 * wa * (1 - a) )
```

#### 2-3. 관련 파생량 계산 루틴 추가
```fortran
! DE 상태방정식
function w_de(a)
  w_de = w0_de + wa_de * (1.0d0 - a)
end function

! DE 밀도 진화
function rho_de(a)
  rho_de = omega_de * a**(-3.0d0*(1.0d0 + w0_de + wa_de)) &
         * exp(-3.0d0 * wa_de * (1.0d0 - a))
end function

! 단열 음속 (c_a² = w - dw/dt / [3H(1+w)])
function ca2_de(a)
  ! dw/da = -wa_de, dw/dt = dw/da * da/dt = dw/da * a * H
  ca2_de = w_de(a) + wa_de * a / (3.0d0 * (1.0d0 + w_de(a)))
end function

! DE Jeans length (c_s ≠ 0 일 때)
function lambda_jeans_de(a)
  lambda_jeans_de = sqrt(cs2_de) / (a * H(a))
end function
```

#### 2-4. H(a)를 사용하는 모든 루틴에 일관 반영
- 시간 ↔ 스케일팩터 변환 (`dt_courant`, `dt_fine` 등)
- kick/drift 연산
- 공동팽창(comoving) 좌표계 관련 함수
- 냉각 함수 내 H(z) 의존 부분
- UV 배경 모델의 적색이동-시간 매핑

#### 2-5. 전처리기
```fortran
#ifdef CPL_DE
  ! CPL 관련 코드
#else
  ! 기존 LCDM 코드 (변경 없음)
#endif
```

#### 2-6. 검증
- `w0_de = -1.0, wa_de = 0.0` 설정 시 기존 ΛCDM 코드와
  bit-wise 동일한 결과 출력 확인

---

## 단계 3-A: c_s² = 0 암흑에너지 섭동 구현

`c_s = 0`이면 DE가 pressureless dust처럼 클러스터링한다.
압력항이 없으므로 실공간에서 직접 풀 수 있다.

### 핵심 방정식 (실공간)

```
∂δ_DE/∂t = -(1+w) * div(v_DE) / a  -  3H * (c_a² - w) * δ_DE

∂v_DE/∂t = -H * v_DE  -  grad(Φ) / a
```

`c_s = 0`이므로 압력 구배항(`c_s² * grad(δ_DE)`)이 없다.

### 수정된 Poisson 방정식

```
∇²Φ = 4πG * a² * ( ρ_m * δ_m  +  ρ_DE(a) * δ_DE )
```

### 구현 상세

#### 3A-1. 격자 변수 추가
AMR 격자에 다음 변수를 추가:
- `delta_de`: DE 밀도 대비 (스칼라)
- `vx_de, vy_de, vz_de`: DE 속도장 (벡터, 3성분)

기존 수력학 변수 배열(`uold`, `unew`)과 유사한 구조로 선언.
또는 별도의 DE 전용 배열 생성.

#### 3A-2. DE 섭동 시간 적분
매 fine 타임스텝마다:

```
Step 1: 현재 Φ로부터 grad(Φ) 계산 (기존 force 루틴 재활용)

Step 2: DE 속도 업데이트 (kick)
        v_DE^(n+1) = v_DE^(n) - dt * [ H * v_DE^(n) + grad(Φ) / a ]

Step 3: DE 밀도 업데이트
        div_v = divergence(v_DE^(n+1))  (유한차분)
        δ_DE^(n+1) = δ_DE^(n) - dt * [ (1+w) * div_v / a
                                       + 3H * (c_a² - w) * δ_DE^(n) ]
```

#### 3A-3. Poisson solver 소스항 수정
기존 multigrid Poisson solver에서 소스항 할당 부분을 찾아 수정:

```fortran
! 기존
rho_source(i,j,k) = rho_matter * delta_matter(i,j,k)

! 수정 (c_s=0, DE 섭동 포함)
#ifdef DE_PERTURB_CS0
rho_source(i,j,k) = rho_matter * delta_matter(i,j,k) &
                   + rho_de(aexp) * delta_de(i,j,k)
#endif
```

`rho_de(aexp)`는 시간에 따라 변하는 DE 배경 밀도.

#### 3A-4. 초기조건
단열(adiabatic) 초기조건 설정:
```
δ_DE(t_init) = (1 + w(a_init)) * δ_m(t_init)
v_DE(t_init) = v_m(t_init)
```

선형 이론에서 c_s=0일 때 DE와 물질의 성장률이 같으므로 이 비례관계가 성립.

#### 3A-5. AMR refinement
- 기존 물질 밀도 기반 refinement criteria에 DE 밀도 기반 조건 추가 (선택사항)
- `delta_de` 필드도 restriction/prolongation 시 처리되어야 함

#### 3A-6. 수치 안정성 주의사항
- `(1 + w)` → 0 근처에서 division by zero 방지
- `w < -1` (phantom regime)에서 `(1+w) < 0` → 부호 반전에 주의
- `c_a²` 계산 시 `(1 + w)` 분모 보호:
  ```fortran
  if (abs(1.0d0 + w) < 1.0d-10) then
    ca2 = w  ! 또는 적절한 limit 값
  else
    ca2 = w + wa_de * aexp / (3.0d0 * (1.0d0 + w))
  endif
  ```

#### 3A-7. 전처리기
```fortran
#ifdef DE_PERTURB_CS0
  ! c_s=0 DE 섭동 코드
#endif
```

---

## 단계 3-B: c_s² ≠ 0 암흑에너지 섭동 구현

`c_s ≠ 0`이면 DE Jeans 스케일이 존재하여 스케일 의존적 클러스터링이 일어난다.
준정적 근사(quasi-static approximation)를 적용하면 DE 섭동이
타원형 방정식(screened Poisson = Helmholtz)으로 환원된다.

이는 RAMSES의 multigrid 구조에 최적이다.

### 핵심 방정식

```
(∇² - 1/λ_J²) * δ_DE = [(1+w) / c_s²] * ∇²Φ
```

여기서 DE Jeans length:
```
λ_J = c_s / (a * H(a))
```

우변의 `∇²Φ`는 Poisson 방정식으로부터:
```
∇²Φ = 4πG * a² * ρ_m * δ_m
```
이므로, 우변은 이미 알고 있는 양이다.

### 구현 상세

#### 3B-1. Helmholtz multigrid solver 생성

기존 multigrid Poisson solver를 복제하여 `helmholtz_fine.f90` (또는 유사한 이름)으로
새 파일을 만들고, relaxation sweep만 수정한다.

**기존 Poisson (Red-Black Gauss-Seidel):**
```fortran
phi_new(i,j,k) = ( phi(i+1,j,k) + phi(i-1,j,k) &
                 + phi(i,j+1,k) + phi(i,j-1,k) &
                 + phi(i,j,k+1) + phi(i,j,k-1) &
                 + dx**2 * source(i,j,k) ) / 6.0d0
```

**Helmholtz 수정:**
```fortran
mass_term = dx**2 / lambda_jeans**2

delta_de_new(i,j,k) = ( delta_de(i+1,j,k) + delta_de(i-1,j,k) &
                       + delta_de(i,j+1,k) + delta_de(i,j-1,k) &
                       + delta_de(i,j,k+1) + delta_de(i,j,k-1) &
                       + dx**2 * source_de(i,j,k) ) &
                     / (6.0d0 + mass_term)
```

분모에 `dx²/λ_J²` 항 하나만 추가된다.
나머지 multigrid 구조(restriction, prolongation, V-cycle)는 완전히 동일하게 유지.

참고: mass term이 있으면 multigrid 수렴이 오히려 개선된다 (스펙트럼 gap 증가).

#### 3B-2. 매 타임스텝 알고리즘

```
Step 1: 물질 밀도장으로 표준 Poisson 풀기
        ∇²Φ⁽⁰⁾ = 4πG * a² * ρ_m * δ_m
        → 기존 Poisson solver 그대로 사용

Step 2: Helmholtz 방정식의 소스항 계산
        source_de(i,j,k) = [(1+w(a)) / cs2_de] * laplacian_phi(i,j,k)
        
        여기서 laplacian_phi는:
        (a) 유한차분으로 직접 계산, 또는
        (b) Step 1의 소스항 자체 (= 4πG a² ρ_m δ_m) 를 재활용
            → 이 방법이 더 정확하고 효율적

Step 3: Helmholtz 방정식 풀기 (새 multigrid solver)
        (∇² - 1/λ_J²) * δ_DE = source_de

Step 4: 수정된 Poisson 방정식 풀기
        ∇²Φ = 4πG * a² * (ρ_m * δ_m + ρ_DE(a) * δ_DE)

Step 5: (선택) 자기일관성 반복
        Step 4의 Φ로 Step 2~4를 반복
        보통 1-2회 반복이면 충분히 수렴
        
Step 6: ∇Φ로 입자 kick + 바리온 수력학 소스항
```

#### 3B-3. AMR 레벨별 최적화 (중요)

각 AMR 레벨의 격자 크기 `dx_level`과 `λ_J`를 비교:

```fortran
! 레벨 ℓ에서
if (dx_level < 0.1d0 * lambda_jeans_de(aexp)) then
  ! 이 레벨에서 DE 섭동이 지수적으로 억제됨
  ! Helmholtz solver 건너뛰고 delta_de = 0 으로 설정
  delta_de(:,:,:) = 0.0d0
else
  ! Helmholtz solver 실행
  call helmholtz_multigrid(...)
endif
```

이 최적화로 세밀한 AMR 레벨(halo 내부 등)에서 불필요한 연산을 생략.

#### 3B-4. λ_J(a) 계산

```fortran
subroutine compute_lambda_jeans(a, lambda_j)
  implicit none
  real(dp), intent(in)  :: a
  real(dp), intent(out) :: lambda_j
  real(dp) :: H_a
  
  call compute_hubble(a, H_a)
  
  ! λ_J = c_s / (a * H)  [comoving 단위]
  lambda_j = sqrt(cs2_de) / (a * H_a)
  
  ! 주의: 코드 내부 단위로 변환 필요
  ! lambda_j = lambda_j / boxlen  (boxlen 단위라면)
end subroutine
```

매 coarse 타임스텝마다 갱신하면 충분 (λ_J는 천천히 변함).

#### 3B-5. 경계조건
- 주기적 경계: 기존 Poisson solver와 동일
- AMR 레벨 간: 기존 prolongation/restriction 연산자 재사용
- 조대 → 세밀 레벨: δ_DE 보간 (기존 인프라 활용)

#### 3B-6. c_s² 값에 따른 자동 분기

```fortran
if (cs2_de < 1.0d-10) then
  ! c_s ≈ 0: Case A 방식 (pressureless fluid)
  call evolve_de_cs0(...)
else if (cs2_de > 0.99d0) then
  ! c_s ≈ 1: 표준 quintessence → δ_DE ≈ 0, 배경만 수정
  ! Helmholtz 풀 필요 없음
  delta_de(:,:,:) = 0.0d0
else
  ! 0 < c_s² < 1: Helmholtz solver 사용
  call evolve_de_helmholtz(...)
endif
```

#### 3B-7. 전처리기
```fortran
#ifdef DE_PERTURB_HELM
  ! c_s≠0 Helmholtz DE 섭동 코드
#endif
```

---

## 단계 4: 새 파일 구조

다음 새 파일들을 생성해줘:

```
ramses/
├── amr/
│   └── ...  (기존, H(a) 수정)
├── poisson/
│   ├── ...  (기존 multigrid, 소스항 수정)
│   └── helmholtz_fine.f90       ← 새 파일: Helmholtz multigrid solver
├── hydro/
│   └── ...  (기존, 최소 수정)
├── pm/
│   └── ...  (기존, kick/drift에 H(a) 반영)
└── dark_energy/                  ← 새 디렉토리
    ├── dark_energy_commons.f90   ← DE 파라미터, 공용 변수
    ├── dark_energy_background.f90 ← w(a), H(a), rho_de(a), lambda_J(a)
    ├── dark_energy_cs0.f90       ← Case A: c_s=0 섭동 진화
    └── dark_energy_helmholtz.f90 ← Case B: c_s≠0 Helmholtz 연동
```

### dark_energy_commons.f90
```fortran
module dark_energy_commons
  use amr_commons, only: dp
  implicit none
  
  ! namelist 파라미터
  real(dp) :: w0_de    = -1.0d0    ! CPL w0
  real(dp) :: wa_de    =  0.0d0    ! CPL wa
  real(dp) :: cs2_de   =  1.0d0    ! DE sound speed squared
  logical  :: de_perturb = .false.  ! DE 섭동 활성화
  
  ! 파생량 (매 스텝 갱신)
  real(dp) :: w_de_current          ! w(a) 현재값
  real(dp) :: rho_de_current        ! ρ_DE(a) 현재값
  real(dp) :: ca2_de_current        ! c_a²(a) 현재값
  real(dp) :: lambda_j_current      ! λ_J(a) 현재값
  
end module
```

Makefile도 새 파일들을 포함하도록 수정해줘.

---

## 단계 5: 검증 테스트

### 5-1. ΛCDM 복원 테스트
```
w0_de = -1.0
wa_de = 0.0
de_perturb = .false.
```
설정으로 기존 ΛCDM 결과와 bit-wise 비교.
H(a), 파워스펙트럼, halo mass function 등.

### 5-2. 배경 우주론 검증
다양한 (w0, wa) 조합에서 H(a), 공동팽창 거리를
CLASS/CAMB 출력과 비교하는 Python 스크립트 작성.

테스트 파라미터:
```
(w0, wa) = (-1.0,  0.0)   # ΛCDM
(w0, wa) = (-0.9,  0.0)   # constant w
(w0, wa) = (-0.9,  0.1)   # thawing
(w0, wa) = (-1.1,  0.2)   # phantom crossing
(w0, wa) = (-0.8, -0.2)   # freezing
```

### 5-3. 선형 섭동 검증
큰 박스(L = 1000 Mpc/h), 낮은 해상도에서 실행.
CLASS/CAMB의 선형 P(k,z)와 RAMSES 결과 비교.

각각 테스트:
```
cs2_de = 0.0    → Case A (DE 클러스터링)
cs2_de = 0.01   → Case B (Helmholtz)
cs2_de = 0.1    → Case B
cs2_de = 1.0    → δ_DE ≈ 0 확인
```

비교 플롯 생성하는 Python 스크립트 포함.

### 5-4. Helmholtz solver 수렴 테스트
- multigrid V-cycle 반복 횟수 대비 잔차(residual) 감소율
- Poisson solver와 비교하여 수렴이 동등하거나 빠른지 확인
- mass term (`dx²/λ_J²`) 크기에 따른 수렴률 변화

### 5-5. namelist 파일 및 실행 스크립트
각 테스트 케이스에 대해:
- RAMSES namelist 파일
- SLURM/PBS 실행 스크립트
- 분석용 Python 스크립트
를 생성해줘.

---

## 단계 6: 성능 최적화 및 병렬화 확인

### 6-1. MPI 병렬화
- `delta_de` 필드의 ghost zone 통신이 올바른지 확인
- Helmholtz solver의 multigrid V-cycle 레벨 간 MPI 통신
  (기존 Poisson solver의 MPI 구조를 동일하게 따르는지)
- Case A의 `v_DE` 필드도 ghost zone 통신 포함

### 6-2. 메모리 오버헤드
- Case A: `delta_de` + `vx_de, vy_de, vz_de` = 격자당 4개 변수 추가
- Case B: `delta_de` = 격자당 1개 변수 추가
- AMR 레벨별 최적화: `dx < 0.1*λ_J`인 레벨에서 DE 배열 할당 생략 가능?

### 6-3. 타이밍
- DE solver 추가에 의한 wall-clock time 증가율 측정
- Helmholtz solve vs Poisson solve 시간 비교
- AMR 레벨 최적화 적용 전/후 비교

프로파일링 결과 요약 스크립트도 작성해줘.

---

## CLAUDE.md (프로젝트 루트에 배치)

아래 내용으로 `CLAUDE.md` 파일을 프로젝트 루트에 생성해줘.
Claude Code가 매 세션마다 자동으로 참조하게 된다.

```markdown
# RAMSES CPL Dark Energy Project

## 목표
RAMSES에 CPL (w0, wa) 동적 암흑에너지 모형 구현
- 배경: H(a) with w(a) = w0 + wa*(1-a)
- c_s=0: DE를 pressureless fluid로 격자 추적
- c_s≠0: 준정적 근사 + Helmholtz multigrid solver

## 코드 규칙
- Fortran 90, 기존 RAMSES 코딩 스타일 준수
- 전처리기: #ifdef CPL_DE (배경), #ifdef DE_PERTURB_CS0 (c_s=0),
            #ifdef DE_PERTURB_HELM (c_s≠0)
- 새 디렉토리: dark_energy/
- 새 파일: dark_energy_commons.f90, dark_energy_background.f90,
          dark_energy_cs0.f90, dark_energy_helmholtz.f90,
          helmholtz_fine.f90
- namelist 파라미터: w0_de, wa_de, cs2_de, de_perturb
- LCDM 호환: w0=-1, wa=0이면 기존 코드와 동일한 결과

## 핵심 방정식
- H²(a) = H0²[Ωm a^(-3) + Ωde a^(-3(1+w0+wa)) exp(-3wa(1-a))]
- w(a) = w0 + wa*(1-a)
- c_a² = w - (a*wa) / [3(1+w)]
- Case A (c_s=0): ∂δDE/∂t = -(1+w)div(vDE)/a - 3H(ca²-w)δDE
- Case B (c_s≠0): (∇² - λJ⁻²)δDE = [(1+w)/cs²]∇²Φ
- λJ = cs/(aH)
- 수정 Poisson: ∇²Φ = 4πGa²(ρm δm + ρDE δDE)

## 수치 안정성
- |1+w| < 1e-10 일 때 보호 코드 필요
- phantom regime (w < -1) 지원
- Helmholtz: dx < 0.1*λJ인 AMR 레벨에서 solver 생략

## 테스트
- w0=-1, wa=0 → LCDM 결과 재현 필수
- CLASS/CAMB P(k) 비교로 선형 영역 검증
- (w0,wa) = (-0.9,0.1), (-1.1,0.2), (-0.8,-0.2) 테스트
- cs2 = 0, 0.01, 0.1, 1.0 테스트
```

---

## 작업 순서 및 git 관리

각 단계 완료 후 반드시 `git commit` 하고 다음 단계로 진행할 것.

| 순서 | 내용 | git branch / commit 메시지 |
|:---:|------|---------------------------|
| 1 | 코드 구조 파악 | (커밋 불필요, 분석만) |
| 2 | 배경 우주론 수정 | `feature/cpl-background` / "Add CPL w0wa background cosmology" |
| 3-A | c_s=0 DE 섭동 | `feature/de-perturb-cs0` / "Add c_s=0 DE perturbation module" |
| 3-B | c_s≠0 Helmholtz | `feature/de-perturb-helm` / "Add Helmholtz multigrid DE solver" |
| 4 | 파일 구조 정리 | 위 브랜치에 포함 |
| 5 | 검증 테스트 | `feature/de-tests` / "Add DE verification tests and scripts" |
| 6 | 최적화 | `feature/de-optimize` / "Optimize DE solver performance" |
