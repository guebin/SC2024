### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 2db4df04-0877-42f7-a2f4-1e7d82da88ea
using LinearAlgebra, PlutoUI,Random

# ╔═╡ 73aa567a-082c-11ef-19ff-8947f8f86252
md"""
# 09wk-2: 실대칭행렬, PD-행렬, SVD의 존재
"""

# ╔═╡ dc9a32c9-366b-453b-a876-d13f8b7d29f3
md"""
## 1. 강의영상
"""

# ╔═╡ eb84030e-e9c1-48c7-8f9e-138b7e4adfee
html"""
<div style="display: flex; justify-content: center;">
<div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
<iframe src=
"
https://youtube.com/embed/playlist?list=PLQqh36zP38-yWPXpE-77ECuOeKaf7U6BC&si=SZZnOiPyZ2U2-upY
"
width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ f9ccb632-23d6-4531-8e3f-b076685baeb6
md"""
## 2. Imports
"""

# ╔═╡ 586aeb22-c4d3-4eb0-8ed7-49ccebb625a1
PlutoUI.TableOfContents()

# ╔═╡ f32588dc-ebeb-4477-8225-045154104adb
md"""
## 3. 실대칭행렬
"""

# ╔═╡ e7c03452-fb03-4907-b689-8e2afd587961
md"""
### A. 스펙트럼 정리 (실수)
"""

# ╔═╡ d0731163-6d63-4abd-a1db-11e6b462559a
md"""
!!! info "이론: 스펙트럼 정리 (Spectral Theorem)" 
	임의의 정사각행렬 ${\bf A}$에 대하여 아래가 성립한다. (사실 정사각행렬 조건을 줄 필요도 없음)

	- ``{\bf A}`` 가 실대칭행렬 $\Leftrightarrow$ ${\bf A}$는 (1) 모든 고유값이 실수이고 (2) 고유벡터가 직교행렬

	여기에서 $\Rightarrow$ 방향을 스펙트럼정리라고 한다. (반대방향은 당연해서..) 
"""

# ╔═╡ a28636dc-dd51-48cf-a3db-23c079459f84
md"""
-- 예제1: 스펙트럼 정리의 역이 성립함을 체크하라.
"""

# ╔═╡ b6d651cf-1131-4d27-86ad-102b8149a31f
md"""
!!! warning "직교행렬의 숨은 조건" 
	``{\bf A}``가 직교행렬이라는 것은 ${\bf A}^{-1}={\bf A}^\top$ 인 경우를 의미한다. 우선 역행렬이 존재해야 하므로 직교행렬이라는 것은 역행렬의 존재성을 암시한다. (즉 full-rank matrix 라는 소리에여) 또한 역행렬이라는 것은 애초에 정사각행렬에 대하여 정의되므로 직교행렬이라는 것은 정사각행렬임을 암시한다. 그리고 직교행렬의 정확한 정의는 
	- 실수인 행렬 ${\bf A}$ 중에서 ${\bf A}^{-1}={\bf A}^\top$ 를 만족하는 행렬
	이므로 직교행렬은 모든 원소가 실수라는 조건을 암시한다. 
"""

# ╔═╡ b0ff6316-aa1e-47b5-a7fc-6a46723b987a
md"""
!!! warning "직교행렬은 왜 실수로 한정했을까?" 
	직교행렬은 (1) 모든 열들이 크기가 1이고 (2) 모든 열들이 서로 수직인 행렬을 상상하기 위해서 고안되었다. 따라서 어떠한 벡터의 (1) 크기와 (2) 각도를 재고 싶은데 실수인 행렬에서 이를 수행하는 적절한 연산은 transpose 이다. 반면에 복소수인 행렬에서 이를 수행하는 적절한 연산은 conjugate transpose 이다. 복소수행렬은 실수행렬을 포함하으로 일반적으로 (1) 모든 열들이 크기가 1이고 (2) 모든 열들이 서로 서로 수직인 행렬의 개념을 잘 구현하는 것은 아래의 조건이다. 

	$${\bf A}^H {\bf A} = {\bf A}{\bf A}^{H} = I$$

	이러한 행렬은 unitary matrix (유니터리 행렬) 이라고 한다. 여기에서 $^H$는 conjugate transpose를 의미한다. 그런데 사실 이런 행렬은 너무 어려워보이고 우리가 별로 관심도 없어서 ${\bf A}$를 실수행렬이라고 가정하고 위의 조건을 다시 써보면 

	$${\bf A}^\top {\bf A} = {\bf A}{\bf A}^{\top} = I$$

	와 같이 되는데, 이러한 ${\bf A}$를 우리는 직교행렬이라고 부른다. 
"""

# ╔═╡ e1c769a5-eced-4276-b9d5-ede43acd5a41
md"""
### B. 스펙트럼 정리의 (복소수)
"""

# ╔═╡ 547f758d-5e11-4da6-87f0-14551d6e8fa8
md"""
!!! info "이론: 스펙트럼 정리 (Spectral Theorem) -- 복소수버전" 
	임의의 정사각행렬 ${\bf A}$에 대하여 아래가 성립한다. (사실 정사각행렬 조건을 줄 필요도 없음)

	- ``{\bf A}`` 가 에르미트행렬 $\Leftrightarrow$ ${\bf A}$는 (1) 모든 고유값이 실수이고 (2) 고유벡터가 유니터리 행렬임. 

	여기에서 $\Rightarrow$ 방향을 스펙트럼정리라고 한다. 
"""

# ╔═╡ ff5714ff-a79e-490e-ae65-bae61999baf5
md"""
!!! tip "스펙트럼 정리의 기억방법" 
	대충 아래와 같이 기억하면 된다. 
	- ``{\bf A}`` 가 실대칭행렬 $\Leftrightarrow$ ${\bf A} = {\bf \Psi}{\bf \Lambda}{\bf \Psi}^\top$ 와 같이 쓸 수 있음. 이때 모든 행렬은 실수임. 
	- ``{\bf A}`` 가 에르미트행렬 $\Leftrightarrow$ ${\bf A} = {\bf \Psi}{\bf \Lambda}{\bf \Psi}^H$ 와 같이 쓸 수 있음. 이때 ${\bf \Lambda}$만 실수임. 
"""

# ╔═╡ 8df6b3c7-6563-475d-88e7-19100f30488e
md"""
-- 예제2. 아래의 행렬들을 관찰하고 스펙트럼 정리가 만족하는지 체크해보자.
"""

# ╔═╡ b8ec122e-c5a8-42c2-b579-d7e12ff62a5a
let 
	# 대칭행렬O, 에르미트행렬O, 스펙트럼정리 만족 (고유값실수 + ΨΛΨ'=A 만족)
	A = [1 3
		 3 2]
	λ,Ψ = eigen(A)
	@show λ # 모든고유값실수O
	@show Ψ'Ψ ≈ Ψ*Ψ' ≈ I # 고유벡터행렬은 직교행렬O
	@show Ψ*Diagonal(λ)*Ψ' ≈ A
end 

# ╔═╡ 192d3686-eb82-4bcb-bb39-54325e3e173c
let 
	# 대칭행렬X, 에르미트행렬O, 스펙트럼정리 만족 (고유값실수 + ΨΛΨ'=A 만족)
	A = [1 -im
		 im  1]
	λ,Ψ = eigen(A)
	@show λ # 모든고유값실수O
	@show Ψ'Ψ ≈ Ψ*Ψ' ≈ I # 고유벡터행렬은 유니터리O
	@show Ψ*Diagonal(λ)*Ψ' ≈ A	
end 

# ╔═╡ 8128fae1-214b-4a4b-822d-ae1f4d670d3a
let 
	# 대칭행렬O, 에르미트행렬X, 스펙트럼정리 만족X (고유값실수X + ΨΛΨ'=A 만족 X)
	A = [im 1
		 1  2]
	λ,Ψ = eigen(A)
	@show λ # 모든고유값실수X 
	@show Ψ'Ψ ≈ Ψ*Ψ' ≈ I # 고유벡터행렬은 유니터리X
	@show Ψ*Diagonal(λ)*Ψ' ≈ A		
end 

# ╔═╡ 4ee89af2-e782-4e25-aba8-40841d5431aa
let 
	# 대칭행렬X, 에르미트행렬X, 스펙트럼정리 만족X (고유값실수X + ΨΛΨ'=A 만족 X)
	A = [1 -1
		 1  2]
	λ,Ψ = eigen(A)
	@show λ # 모든고유값실수X
	@show Ψ'Ψ ≈ Ψ*Ψ' ≈ I # 고유벡터행렬은 유니터리X
	@show Ψ*Diagonal(λ)*Ψ' ≈ A		
end 

# ╔═╡ 9a1f6b57-d0a5-417b-9f25-eeeb920cfb3a
let 
	# (그래프) 퓨리에변환 연구할때 잘나오는 행렬.. 
	# 대칭행렬X, 에르미트행렬X, 스펙트럼정리 만족X (고유값실수X + ΨΛΨ'=A 만족 X)
	A = [0 1 0 
	     0 0 1 
	     1 0 0]
	λ,Ψ = eigen(A)
	@show λ # 모든고유값실수X
	@show Ψ'Ψ ≈ Ψ*Ψ' ≈ I # 고유벡터행렬은 유니터리O
	@show Ψ*Diagonal(λ)*Ψ' ≈ A			
end 

# ╔═╡ 0d3ae451-7512-4016-8394-843091a71767
md"""
### C. 생각해볼 점들 
"""

# ╔═╡ 62c5757a-6807-4c73-9b7f-88128b0d4370
md"""
-- 예제1. 실대칭행렬일때 "고유값분해 = 특이값분해" 가 성립하지 않는 이유를 설명해보라. 
"""

# ╔═╡ 701374e0-21be-4db6-aea4-77b8166f3909
# 고유값이 음수 나올수 있으니까

# ╔═╡ 7b6c2413-79c7-491c-b872-d15fa9944685
md"""
-- 예제2. 실대칭행렬 ${\bf A}$의 고유값중 0이 하나도 없다면, ${\bf A}$의 역행렬이 존재함을 보여라. (사실 실대칭행렬 아니고 대각화가능행렬이기만 해도 가능함)
"""

# ╔═╡ 470c3022-480f-4291-9b9b-0ff8f6efe694
let 
	A = [1 3
		 3 2]
	λ,Ψ = eigen(A)
	@show λ # 고유값에 0은 없음 
	λ1,λ2 = λ
	Ainv = Ψ * Diagonal([1/λ1, 1/λ2]) * Ψ'
	@show A*Ainv ≈ Ainv*A ≈ I
end 

# ╔═╡ d480d5f3-53ae-4fb6-877c-60206626f80f
md"""
## 4. PD행렬
"""

# ╔═╡ 17fc6337-1ad4-40bf-8887-9930ec76a21e
md"""
### A. PD/PSD, ND/NSD matrix의 정의
"""

# ╔═╡ 06f5d756-8656-431d-834d-b6f591fe4967
md"""
!!! info "정의 : PD/PSD, ND/NSD" 
	임의의 non-zero ${\bf y} \in \mathbb{R}^n$ (non-zero ${\bf y}\in \mathbb{C}^n$ ) 에 대하여 아래를 만족하는 실대칭행렬 (혹은 에르미트행렬) ${\bf A}_{n \times n}$를 positive definite matrix 라고 한다. 

	${\bf y}^\top {\bf A} {\bf y} > 0 \quad \quad ({\text or}\quad {\bf y}^H {\bf A} {\bf y} > 0)$

	이외에도 $\geq, < ,\leq$ 에 따라 postive semi-definite, negative definite, negative semi-definite matrix를 정의할 수 있다. 
"""

# ╔═╡ 1ec3be98-d795-425e-8522-57565b0d116f
md"""
-- 예제1: PD-행렬인 경우
"""

# ╔═╡ dae2fe6b-780a-4875-9164-c6de05b9fcca
let 
	A = [1.0 0.5
		 0.5 1.0]
	y = [15,-22]
	y'A*y	
end 

# ╔═╡ d588a981-351c-44dc-ac96-829e953d61a3
md"""
-- 예제2: 실대칭행렬(에르미트)이지만 PD-행렬이 아닌경우
"""

# ╔═╡ b24059c2-1205-49a0-8bfb-c0b4537c2da4
let 
	A = [1 3
		 3 2]
	y = [15,-22]
	y'A*y
end 

# ╔═╡ 1d247622-5ab9-4207-9b49-35e5c0eae9e4
md"""
### B. PD-행렬의 성질
"""

# ╔═╡ 21588176-bc55-4d79-82b2-6192e3871ffc
md"""
!!! info "정리: PD행렬 = 스펨트럼정리를 만족 + 모든고유값이 양수" 
	임의의 정사각행렬 ${\bf A}$에 대하여 아래가 성립한다. (사실 정사각행렬 조건을 줄 필요도 없음)
	- ``{\bf A}`` 가 PD-행렬 $\Leftrightarrow$ (1) 모든 고유값이 양수 (2) 고유벡터행렬이 직교행렬 (혹은 유니터리 행렬)

	``{\bf A}`` 가 PSD, ND, NSD 인 행렬에 대해서도 비슷한 정리가 성립한다. 
"""

# ╔═╡ cd03df6c-5473-4648-a51c-09b2438078b9
md"""
-- 예제1: PD-행렬인 경우
"""

# ╔═╡ f2ad8b69-812c-481c-a590-76dd8d81f774
let 
	# A는 실대칭행렬 & 모든고유값이 양수 -> PD-행렬
	A = [1.0 0.5
		 0.5 1.0]
	eigen(A)
end 

# ╔═╡ f923156f-709f-40af-8c45-b80c7e15383c
md"""
-- 예제2: 모든 non-zero real vector ${\bf y}$에 대하여 ${\bf y}^\top {\bf A}{\bf y}>0$ 이지만, PD-행렬이 아닌 경우 (이유? 실대칭행렬 X)
"""

# ╔═╡ f5177758-e2c0-42cc-a11e-f7f3b2935378
y1 = @bind y1 Slider(-50:0.1:50,show_value=true)

# ╔═╡ 6807722e-1124-4fd1-8e5f-22771f534d56
y2 = @bind y2 Slider(-50:0.1:50,show_value=true)

# ╔═╡ b234b344-7577-43a1-a4e5-bdaf28620d4d
let
	A = [1 -1 
		 1  1]
	y = [-2.5, -10.5]
	@show y'A*y # 항상양수
	λ,Ψ = eigen(A)
	@show λ # 모든고유값실수X
	@show Ψ'Ψ ≈ Ψ*Ψ' ≈ I # 고유벡터행렬은 유니터리
	@show Ψ*Diagonal(λ)*Ψ' ≈ A # 이거 가능
end 

# ╔═╡ 1fa71490-f3bb-433f-a606-99b1217655d8
md"""
-- 예제3: 모든 non-zero real vector ${\bf y}$에 대하여 ${\bf y}^\top {\bf A}{\bf y}>0$ 이고 모든 고유값도 양수이지만, PD-행렬이 아닌 경우 (이유? 실대칭행렬 X)
"""

# ╔═╡ b275b633-ee83-4473-bb00-e05552697a59
let 	
	A = [1 1
		 0 2]
	y = [-2.5, -10.5]
	@show y'A*y # 항상양수	
	λ,Ψ = eigen(A)
	@show λ # 모든고유값실수O + 심지어 양수
	@show Ψ'Ψ ≈ Ψ*Ψ' ≈ I # 고유벡터행렬은 직교행렬 X
	@show Ψ*Diagonal(λ)*Ψ' ≈ A # 이거 불가능
	Ψ
end 

# ╔═╡ f6e131e4-d5a6-488c-bb71-78aa9ab09425
md"""
### C. 생각해 볼 점
"""

# ╔═╡ f25cb325-a290-4156-b978-c4de13840b20
md"""
-- 예제1: 임의의 실수행렬 ${\bf X}_{n\times p}$에 대하여 ${\bf X}^\top{\bf X}$ 와 ${\bf X}{\bf X}^\top$는 PD-행렬임을 보여라.
"""

# ╔═╡ 865e099e-8c41-45a3-b43e-c8bc820ad979
# X'X 
## 대칭행렬임을 보이는건 쉬움: (X'X)' = X'X
## ∀y: y'Ay >0 임을 보이면 되는데 y'X'X'y = ỹ'ỹ > 0  

# XX' 
## 생략.. 

# ╔═╡ 42a677cc-06b4-473c-be9d-00ecee71448e
md"""
-- 예제2: 임의의 행렬 ${\bf A}$가 PD-행렬이라면 아래를 만족하는 행렬들이 존재함을 보여라. (찾으면 됩니당)

- ``{\bf A}^{1/2}{\bf A}^{1/2}={\bf A}``
- ``{\bf A}^{-1/2}{\bf A}^{-1/2}={\bf A}^{-1}``
"""

# ╔═╡ 18994f7f-c5ce-4a3d-8131-338d027022d1
let 
	# 이런논리로 찾으면 된다. 
	A = [1.0 0.5
		 0.5 1.0]
	λ,Ψ = eigen(A)
	λ1,λ2 = λ 
	B = Ψ * Diagonal([√λ1, √λ2]) * Ψ' # B는 √A 같은거...
	C = Ψ * Diagonal([1/√λ1, 1/√λ2]) * Ψ' # C는 (√A)^(-1) 같은거..
	@show B*B ≈ A
	@show C*C ≈ inv(A)
end 

# ╔═╡ 1186d3d4-c5cc-4eeb-8d44-eae2aaabd43a
md"""
## 5. 정리
"""

# ╔═╡ 18c16498-017d-4053-81d3-fb1aa9dfefa8
md"""
### A. 여러가지 행렬들...
"""

# ╔═╡ bc77d15c-cba2-4cde-8d77-eb009df05731
let 
	A = [1 2
		 0 -3]
	y = [1,-1]
	y'A*y

end 

# ╔═╡ d2598cad-acdd-4cf6-aca3-44ee2f846c36
md"""
| 예시 | (1)``{\bf A}{\bf \Psi}= {\bf \Psi}{\bf \Lambda}`` | (2)``{\bf A}={\bf \Psi}{\bf \Lambda}{\bf \Psi}^{-1}`` | (3)``{\bf A}={\bf \Psi}{\bf \Lambda}{\bf \Psi}^H`` | (4)``\forall \lambda_i \in \mathbb{R}`` | (5)``\forall\lambda_i  >0`` | (6)``{\bf y}^H {\bf A}{\bf y}>0`` |
|:---------:|:---------:|:---------:|:---------:|:--------:|:---------:|:---------:|
|``\begin{bmatrix} i & i \\ 0 & i \end{bmatrix}`` |  O   | X  | X   | X  | NA | X |
|``\begin{bmatrix} -1 & -1 \\ 0 & -1 \end{bmatrix}`` |  O   | X  | X   | O | X | X |
|``\begin{bmatrix} 1 & 2 \\ 0 & -3 \end{bmatrix}`` |  O   | O  | X   | O  | X | X |
|``\begin{bmatrix} 0 & 1 & 0 \\ 0 & 0 & 1 \\ 1 & 0 & 0 \end{bmatrix}`` |  O   | O  | O |  X  | NA | X |
|``\begin{bmatrix} 1 & 3 \\ 3 & 2 \end{bmatrix}`` |  O   | O  | O  | O  | X | X |
|``\begin{bmatrix} 1 & i \\ -i & -3 \end{bmatrix}`` |  O   | O  | O  | O  | X | X |
|``\begin{bmatrix} 1 & 0.5 \\ 0.5 & 1 \end{bmatrix}`` |  O   | O  | O  | O  | O | O|
|``\begin{bmatrix} 2 & i \\ -i & 2 \end{bmatrix}`` |  O   | O  | O  | O  | O | O|
|``\begin{bmatrix} 1 & -1 \\ 1 & 1 \end{bmatrix}`` |  O   | O  | O  | X  | NA | O|
|``\begin{bmatrix} 1 & 1 \\ 0 & 2 \end{bmatrix}`` |  O   | O  | X  | O  | O | O|
"""

# ╔═╡ 0b153dcc-cb5e-4651-832b-697353ec8dc1
md"""
-- 잔소리1
- (1)은 정사각행렬이기만 하면 성립한다. 
- (2)는 고유벡터행렬이 full-rank 일 경우 성립한다. 
- 만약 모든 고유값이 서로 다르면 (2)는 무조건 성립한다. 
- (3)이 성립하면 자동으로 (2)가 성립된다. 
- (1)-(3)이 성립해도 (4)가 성립하지 않을수 있다. (GFT)
- 에르미트행렬이면 (1)-(4)가 무조건 성립한다. (스펙트럼 정리) 
- (3)-(4)가 성립하면 (혹은 (1)-(4) 가 성립하면) 에르미트행렬이다. 
- 에르미트행렬이면 (5),(6)은 동시에 성립하거나 동시에 성립하지 않는다. 즉 에르미트행렬이라는 전제를 주면 (5)와 (6)은 동치이다. 
- (6)이 성립한다고 고유값이 항상 양수인 것은 아니다. (마지막 두번째예시)
- (5),(6)이 동시에 성립하지만 고유벡터가 직교하지 않을 수 있다. (마지막예시)
- **(1)-(6)의 조건에서 "고유값분해 = 특이값분해" 가 성립한다.** <-- 중요해요
"""

# ╔═╡ c86ab7ca-7c75-4016-872f-e7447d64a277
md"""
-- 잔소리2
- (1) 을 만족하는 행렬은 정사각행렬임. 
- (1)-(2) 를 만족하는 행렬을 대각화가능행렬이라고 함. 
- (1)-(4) 를 만족하는 행렬은 (애매하지만) 직교대각화가능행렬 이라고 하는 경우가 많음. 즉 스펙트럼 정리를 만족하는 행렬을 직교대각화가능행렬이라고 표현함. 
- (1)-(6) 을 만족하는 행렬은 PD-행렬이라고 표현함. 
"""

# ╔═╡ 1de9979f-6f2d-4f40-8038-c52172a8ad15
md"""
!!! warning "공부할 때 조심할 점"
	교재나 문헌들이 은근히 엄밀하지 않게 쓰는 경우가 많습니다. 그래서 따져봐야해요.. 
	- <https://dept.math.lsa.umich.edu/~kesmith/SpectralTheorem.pdf>
	- <https://math.emory.edu/~lchen41/teaching/2020_Fall/Section_8-2.pdf>
	- <https://en.wikipedia.org/wiki/Spectral_theorem>
"""

# ╔═╡ 54dc8a51-192d-45ae-bb73-a0afad4aad28
md"""
### B. 행렬의 숨겨진 조건
"""

# ╔═╡ 3cdf9793-c896-432b-ac17-ba49965e57f1
md"""
!!! warning "대칭행렬, 에르미트행렬, 실대칭행렬"
	아래는 교재/문헌을 읽으면서 눈치껏 알아야하는 내용들을 정리한 것 이다. 
	1. ``{\bf A}``가 대칭행렬이라는 의미는 정사각행렬임을 암시한다. 아무말없이 ``{\bf A}``가 대칭행렬이라고만 나와도 알아서 정사각행렬이라 생각해야 한다.
	2. 당연히 에르미트행렬과 실대칭행렬도 정사각행렬을 전제하고 읽어야 한다. 
	3. 책을 읽다가 에르미트행렬이라는 단어가 나오면 원소가 복소수라고 생각하는 편이 낫다. (실제로 에르미트행렬의 정의가 그렇지 않다고 하더라도..) 만약 책을 쓰는 사람이 "모든 원소가 실수인 에르미트 행렬"을 생각하고 싶었다면, 그냥 "실대칭행렬"이라고 지칭했을 것이다. 
	4. 블로그등에서 ``{\bf A}``가 대칭행렬이라고 하는 경우, 실제로는 문맥상 **"실대칭행렬"을 지칭하는 경우가 많다.** 따져보면 실대칭행렬이 아니지만 대칭행렬이라는 것은 에르미트 행렬이 아니란 소리인데, (즉 대칭행렬O / 실대칭행렬X / 에르미트행렬X 인 경우인데) 이러한 행렬은 특성을 연구할 것이 별로 없고 블로그에 포스팅 할것도 별로 없다. 
	행렬별로 정리하면 아래와 같다. 
	1. 실대칭행렬 $\to$ 정사각행렬 
	1. 대칭행렬 $\to$ 정사각행렬 + 원소는 아마도 실수
	2. 에르미트행렬 $\to$ 정사각행렬 + 원소는 아마도 복소수
"""

# ╔═╡ 185376c2-5ad5-4e31-a1d9-04029024d5c3
md"""
!!! warning "직교행렬, 유니터리행렬"
	아래의 내용들도 눈치껏 이해해야한다. 
	1. 직교행렬 $\to$ 모든 원소 실수 + 역행렬존재 + 정사각행렬 
	2. 유니터리행렬 $\to$ 원소는 아마도 복소수 + 역행렬존재 + 정사각행렬 
"""

# ╔═╡ e187b01e-45a2-4457-948b-3a4bba05e9fb
md"""
!!! warning "PD-행렬, PSD-행렬"
	아래의 내용들도 눈치껏 이해해야한다. 
	1. PD-행렬 / PSD-행렬 $\to$ 실대칭행렬 or 에르미트행렬 + 정사각행렬
	2. 가끔씩 PD-행렬인데 대칭을 전제하지 않는 경우도 있음. (무시할 것)
"""

# ╔═╡ 73512dff-60dc-47b2-bac6-7d2faebb3656
md"""
## 6. SVD의 존재성 ($\star$)
"""

# ╔═╡ 2d15f25e-bfe1-4787-830a-f0633e98828b
md"""
이말하려고 먼길 돌아왔습니다. 임의의 행렬 ${\bf X}$에 대하여 SVD는 무조건 존재합니다. 
"""

# ╔═╡ 26059dc2-ecb9-4f8e-845d-a4e0cee2bbc9
md"""
!!! info "정리: 특이값분해 (SVD)" 
	임의의 행렬 ``{\bf X}_{n\times p}`` 에 대하여 특이값 분해는 항상 존재한다. 
"""

# ╔═╡ 75776417-3359-4e3a-8f28-76d7b0a316a6
md"""
-- 예제1: ${\bf V}$의 열들은 ${\bf X}^\top{\bf X}$의 고유벡터임을 체크하라. 이때 대응하는 고유값은 ${\bf D}^2$임을 체크하라. 
"""

# ╔═╡ fd0dc0f6-862c-4dda-a76b-56cfad59f56c
let 
	X = randn(100,5)
	U,d,V = svd(X)
	(X'X)* V ≈ V * Diagonal(d.^2)
end 

# ╔═╡ 68f8a231-03bf-47cb-a86e-3b4457a07259
md"""
-- 예제2: ${\bf U}$의 열들은 ${\bf X}{\bf X}^\top$의 고유벡터임을 체크하라. 이때 대응하는 고유값은 ${\bf D}^2$임을 체크하라. 
"""

# ╔═╡ 10205523-f237-4d41-92ee-eb2d7803cde5
let 
	X = randn(100,5)
	U,d,V = svd(X)
	(X*X') * U ≈ U * Diagonal(d.^2) 
end 

# ╔═╡ 99428e3c-1e8f-408e-bcb1-0b266314becf
md"""
## 7. 특이값분해 = 고유값분해
"""

# ╔═╡ 7f04c611-6835-478d-a4a0-244f63837f2f
md"""
!!! info "정리: 특이값분해 = 고유값분해" 
	행렬 ``{\bf A}`` 가 PSD-행렬이라면 특이값분해와 고유값분해가 일치하도록 만들 수 있다. 즉 아래를 만족하는 직교행렬 ``{\bf U}={\bf V}={\bf \Psi}`` 와 대각원소가 비음인 대각행렬 ``{\bf D}={\bf \Lambda}`` 가 존재한다.

	$${\bf A} = {\bf U}{\bf D}{\bf V}^\top = {\bf \Psi}{\bf \Lambda}{\bf \Psi}^\top$$

	이때 ``{\bf U}={\bf V}={\bf \Psi}``는 ${\bf A}$의 고유벡터행렬이면서 동시에 ${\bf A}^\top{\bf A}$, ${\bf A}{\bf A}^\top$ 의 고유벡터행렬이 된다. 그리고 ${\bf D}={\bf \Lambda}$는 ${\bf A}$의 고유값행렬이면서 동시에 ${\bf A}^\top{\bf A}$, ${\bf A}{\bf A}^\top$ 의 고유값행렬의 $\sqrt{}$ 가 된다.
"""

# ╔═╡ f1d33402-0a82-4615-8df6-feb87ba64944
let 
	Random.seed!(1)
	X = randn(100,2)
	U,d,V = svd(X'X)
	λ,Ψ = eigen(X'X,sortby= - )
	Ψ1,Ψ2 = eachcol(Ψ)
	Ψ = [Ψ1 -Ψ2] # U,V 와 일치하도록 조정 
	@show U ≈ V ≈ Ψ
	@show λ ≈ d
	[U V Ψ]
end

# ╔═╡ 1ecc6b50-8a3c-4af9-bb12-5cc6c1d2be98


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.3"
manifest_format = "2.0"
project_hash = "7f53ddf024047cd4b795dc4ba6acfc656a7c2793"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─73aa567a-082c-11ef-19ff-8947f8f86252
# ╟─dc9a32c9-366b-453b-a876-d13f8b7d29f3
# ╟─eb84030e-e9c1-48c7-8f9e-138b7e4adfee
# ╟─f9ccb632-23d6-4531-8e3f-b076685baeb6
# ╠═2db4df04-0877-42f7-a2f4-1e7d82da88ea
# ╠═586aeb22-c4d3-4eb0-8ed7-49ccebb625a1
# ╟─f32588dc-ebeb-4477-8225-045154104adb
# ╟─e7c03452-fb03-4907-b689-8e2afd587961
# ╟─d0731163-6d63-4abd-a1db-11e6b462559a
# ╟─a28636dc-dd51-48cf-a3db-23c079459f84
# ╟─b6d651cf-1131-4d27-86ad-102b8149a31f
# ╟─b0ff6316-aa1e-47b5-a7fc-6a46723b987a
# ╟─e1c769a5-eced-4276-b9d5-ede43acd5a41
# ╟─547f758d-5e11-4da6-87f0-14551d6e8fa8
# ╟─ff5714ff-a79e-490e-ae65-bae61999baf5
# ╟─8df6b3c7-6563-475d-88e7-19100f30488e
# ╠═b8ec122e-c5a8-42c2-b579-d7e12ff62a5a
# ╠═192d3686-eb82-4bcb-bb39-54325e3e173c
# ╠═8128fae1-214b-4a4b-822d-ae1f4d670d3a
# ╠═4ee89af2-e782-4e25-aba8-40841d5431aa
# ╠═9a1f6b57-d0a5-417b-9f25-eeeb920cfb3a
# ╟─0d3ae451-7512-4016-8394-843091a71767
# ╟─62c5757a-6807-4c73-9b7f-88128b0d4370
# ╠═701374e0-21be-4db6-aea4-77b8166f3909
# ╟─7b6c2413-79c7-491c-b872-d15fa9944685
# ╠═470c3022-480f-4291-9b9b-0ff8f6efe694
# ╟─d480d5f3-53ae-4fb6-877c-60206626f80f
# ╟─17fc6337-1ad4-40bf-8887-9930ec76a21e
# ╟─06f5d756-8656-431d-834d-b6f591fe4967
# ╟─1ec3be98-d795-425e-8522-57565b0d116f
# ╠═dae2fe6b-780a-4875-9164-c6de05b9fcca
# ╟─d588a981-351c-44dc-ac96-829e953d61a3
# ╠═b24059c2-1205-49a0-8bfb-c0b4537c2da4
# ╟─1d247622-5ab9-4207-9b49-35e5c0eae9e4
# ╟─21588176-bc55-4d79-82b2-6192e3871ffc
# ╟─cd03df6c-5473-4648-a51c-09b2438078b9
# ╠═f2ad8b69-812c-481c-a590-76dd8d81f774
# ╟─f923156f-709f-40af-8c45-b80c7e15383c
# ╠═f5177758-e2c0-42cc-a11e-f7f3b2935378
# ╠═6807722e-1124-4fd1-8e5f-22771f534d56
# ╠═b234b344-7577-43a1-a4e5-bdaf28620d4d
# ╟─1fa71490-f3bb-433f-a606-99b1217655d8
# ╠═b275b633-ee83-4473-bb00-e05552697a59
# ╟─f6e131e4-d5a6-488c-bb71-78aa9ab09425
# ╟─f25cb325-a290-4156-b978-c4de13840b20
# ╠═865e099e-8c41-45a3-b43e-c8bc820ad979
# ╟─42a677cc-06b4-473c-be9d-00ecee71448e
# ╠═18994f7f-c5ce-4a3d-8131-338d027022d1
# ╟─1186d3d4-c5cc-4eeb-8d44-eae2aaabd43a
# ╟─18c16498-017d-4053-81d3-fb1aa9dfefa8
# ╠═bc77d15c-cba2-4cde-8d77-eb009df05731
# ╟─d2598cad-acdd-4cf6-aca3-44ee2f846c36
# ╟─0b153dcc-cb5e-4651-832b-697353ec8dc1
# ╟─c86ab7ca-7c75-4016-872f-e7447d64a277
# ╟─1de9979f-6f2d-4f40-8038-c52172a8ad15
# ╟─54dc8a51-192d-45ae-bb73-a0afad4aad28
# ╟─3cdf9793-c896-432b-ac17-ba49965e57f1
# ╟─185376c2-5ad5-4e31-a1d9-04029024d5c3
# ╟─e187b01e-45a2-4457-948b-3a4bba05e9fb
# ╟─73512dff-60dc-47b2-bac6-7d2faebb3656
# ╟─2d15f25e-bfe1-4787-830a-f0633e98828b
# ╟─26059dc2-ecb9-4f8e-845d-a4e0cee2bbc9
# ╟─75776417-3359-4e3a-8f28-76d7b0a316a6
# ╠═fd0dc0f6-862c-4dda-a76b-56cfad59f56c
# ╟─68f8a231-03bf-47cb-a86e-3b4457a07259
# ╠═10205523-f237-4d41-92ee-eb2d7803cde5
# ╟─99428e3c-1e8f-408e-bcb1-0b266314becf
# ╟─7f04c611-6835-478d-a4a0-244f63837f2f
# ╠═f1d33402-0a82-4615-8df6-feb87ba64944
# ╠═1ecc6b50-8a3c-4af9-bb12-5cc6c1d2be98
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
