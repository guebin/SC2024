### A Pluto.jl notebook ###
# v0.19.4

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

# ╔═╡ 15442a71-00fe-48ed-ae51-def2f2d8876e
using LinearAlgebra, PlutoUI, Plots

# ╔═╡ e08c8300-db18-11ec-3f60-19c9638fa2d2
md"""
# 5월24일
"""

# ╔═╡ 90db85a9-9864-4e63-98ca-7edb77f1932f
html"""
<div style="display: flex; justify-content: center;">
<div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
<iframe src=
"
https://www.youtube.com/embed/playlist?list=PLQqh36zP38-zjEXJRTGfKk4kWO5T4IlYy
"
width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ c84fb7f0-e572-46ab-b0a3-7612c1870815
md"""
## usings
"""

# ╔═╡ 8308a4c4-46b0-4407-aba5-d8b174c2f152
Plots.plotly()

# ╔═╡ 8b5ce365-94e4-4875-ac48-d7d6f4b7409c
PlutoUI.TableOfContents()

# ╔═╡ 9d9020be-5c18-48f6-bf87-d3623f568673
md"""
## 지난시간 summary
"""

# ╔═╡ 05c1b786-159f-4ebd-9748-65c5194133f6
md"""
### 몇 가지 용어의 약속 (책에 없는 용어도 있음)
"""

# ╔═╡ 65396179-2be2-4ba2-84a6-38db6a4335c5
md"""
(1) ${\bf X}$의 svd 분해꼴로 얻어지는 결과를 각각 ${\bf X}$의 u-matrix, d-matrix, v-matrix 라고 부르자. 
- 참고로 ${\bf X}$의 u-matrix, v-matrix 유일하지 않음. 
- u-matrix와 v-matrix의 column들의 순서를 임의로 바꿀수 있으므로 ${\bf X}$의 d-matrix 역시 유일하진 않음. 그렇지만 순서를 정렬한다면 유일함. 
"""

# ╔═╡ 318c6975-764e-47e5-abe7-0ba7b9883d20
md"""
(코드예시)
"""

# ╔═╡ a411e2b7-b659-44d1-8146-16d15f59e0e4
let 
	X = [1 2 ; 3 4]
	U,d,V = svd(X)
	U # u-matrix of X 
	Diagonal(d) # d-matrix of X 
	V # v-matrix of X
end 

# ╔═╡ 867a6ec4-fb9c-4001-910f-9080ecfcddbc
md"""
(2) 정사각행렬 ${\bf A}$가 고유값분해가 가능할때 ${\bf A}$의 서로독립인 고유벡터를 column-wise하게 합친 행렬을 ${\bf A}$의 **고유벡터행렬(eigenvector matrix)**이라고 하자. 그리고 고유벡터행렬의 each columns에 대응하는 고유값을 순서대로 나열하여 대각선의 원소에 넣은 대각행렬을 **고유값행렬(eigenvalue matrix)**라고 하자. 
- 참고로 ${\bf X}$의 고유벡터행렬은 유일하지 않음. 예를들어 매트릭스 ${\bf \Psi}$가 ${\bf X}$의 고유벡터행렬이면 $-{\bf \Psi}$ 역시 ${\bf X}$의 고유벡터행렬임. 
- 고유벡터의 순서를 바꿀수 있기에 고유값행렬도 유일하지는 않음. 하지만 고유값을 순서대로 정리한다면 유일함.

"""

# ╔═╡ 67b61ca4-a06b-4f1d-8cab-0d4ba12dfae6
md"""
(코드예시)
"""

# ╔═╡ 11d55eb6-67af-4146-abe0-c34dba602e71
let 
	A = [1 2; 2 1]
	λ, Ψ = eigen(A)
	Ψ # eigenvector-matrix of A 
	Diagonal(λ) # eigenvalue-matrix of A 
end

# ╔═╡ 6f9a1a8e-717b-4953-804a-e8647743fbd4
md"""
### 지난시간 요약
"""

# ╔═╡ 598d1691-16cf-4803-a5af-d87abec1a02e
md"""
#### 지난시간까지의 내용 요약
"""

# ╔═╡ fa178a75-8119-4d2c-9f63-94420591e7a5
md"""
(1) 임의의 매트릭스 ${\bf X}_{n\times p}$는 항상 SVD분해가 가능하다. 
"""

# ╔═╡ c1e1aded-6005-4b01-bab1-8f75a3811dc4
md"""
(2) ${\bf X}_{n \times p}$가 이미지일때: SVD를 이용하면 ${\bf X}_{n\times p} \approx {\bf \hat{X}}_{n\times p}$ 인 ${\bf \hat{X}}_{n\times p}$를 구할 수 있다. (압축하니)
"""

# ╔═╡ 752ebf31-3ec1-4033-a3d1-0bddd8ac3232
md"""
(3) ${\bf X}_{n \times p}$가 데이터프레임일때: SVD를 이용하면 ${\bf X}_{n \times p}$보다 rank가 작은 매트릭스 ${\bf Z}_{n \times q <p}$ 를 찾을 수 있다. ${\bf Z}_{n\times q}$는 적당한 변환 ${\bf B}_{q\times p }$에 의하여 ${\bf Z}_{n\times q}{\bf B}_{q\times p}={\bf \hat{X}}_{n\times p} \approx {\bf X}_{n\times p}$ 이 된다는 특징이 있다. 
"""

# ╔═╡ 3a85c442-10cd-46ac-943f-6912018b008c
md"""
#### ${\bf Z}$와 ${\bf B}$에 대한 암기사항
"""

# ╔═╡ 330e2f3a-4d66-417f-aaf5-b92ba0a79c55
md"""
(1) ${\bf X}={\bf U}{\bf D}{\bf V}^\top$ 일때 ${\bf Z}={\bf \tilde U}{\bf \tilde D}={\bf X}{\bf \tilde V}$ 와 같이 구할 수 있다. 
- 이때 ${\bf \tilde U}$, ${\bf \tilde D}$, ${\bf \tilde V}$의 정의는 지난강의노트 참고
"""

# ╔═╡ 713775f6-ccce-4b7a-bda3-e828aebcd6a1
md"""
(2) ${\bf X}$의 v-martrix를 구하는 방법이 다양하다. 
- 방법1: ${\bf X}$의 svd를 이용하여 직접구한다. 
- 방법2: ${\bf X}'{\bf X}$에 svd를 수행한다. 그러면 신기하게도 ${\bf X}'{\bf X}$의 v-matrix와 u-matrix는 항상 서로 값이 같게 나온다. (왜?) 그리고 더 신기하게도 이때의 ${\bf X}'{\bf X}$의 v-matrixs를 (혹은 u-matrix를)  ${\bf X}$의 v-matrix라고 주장할 수 있다. (왜?)
- 방법3: ${\bf X}'{\bf X}$의 고유벡터행렬을 구하면 그 행렬을 ${\bf X}$의 v-matrix라고 주장할 수 있다. (왜?)
"""

# ╔═╡ 006740e0-f12c-4312-b181-30d5e6f82ec3
md"""
## Eigenvalue Decomposition
"""

# ╔═╡ cba26403-d491-4eae-be95-7a11984c2b0e
md"""
### 고유값과 고유벡터의 정의
"""

# ╔═╡ fe5876ba-30f9-470d-bd42-4c188fa580f3
md"""
`-` 임의의 정사각행렬 ${\bf A}_{n\times n}$에 대하여 어떠한 벡터 ${\psi}_{n\times 1}\neq 0$ 가 적당한한 값 $\lambda$에 대하여 

$${\bf A}{\psi} = \lambda \psi$$

를 만족하면 $\psi$를 $\lambda$의 고유벡터라고 하고 $\lambda$는 $\psi$에 대응하는 고유값이라고 한다. 
- note: 0-벡터는 고유벡터로 인정하지 않음 $\to$ 고유값을 찾는방법? $\det({\bf A}-\lambda {\bf I})=0$을 만족하는 $\lambda$를 풀면된다. 
"""

# ╔═╡ bca97f96-16d9-4c2d-8d6b-6bd5c7400706
md"""
(예제1) 첫번째 고유값
"""

# ╔═╡ ee444ee9-9a95-4b8f-999d-60d1c16df56f
let 
	A = [1 2 ; 2 1]
	I = [1 0 ; 0 1]
	λ = -1 
	det(A-λ*I) 
end

# ╔═╡ 4a687148-264d-46b4-a1e4-7d3e94288e57
md"""
-  $\lambda = -1$일때 $\det ({\bf A}-\lambda {\bf I})=0$ 이므로 $\lambda =-1$은 ${\bf A}$의 고유값이다.
"""

# ╔═╡ 17bce861-3ef6-4050-bc79-67da0ea689b8
md"""
(예제2) 두번째 고유값
"""

# ╔═╡ 6349ab90-43fd-4a19-87d2-aa2483d3faea
let 
	A = [1 2 ; 2 1]
	I = [1 0 ; 0 1]
	λ = 3
	det(A-λ*I) 
end

# ╔═╡ d7d62756-2dab-4a0f-817a-e1acd83a3b74
md"""
-  $\lambda = 3$일때 $\det ({\bf A}-\lambda {\bf I})=0$ 이므로 $\lambda =3$은 ${\bf A}$의 고유값이다.
"""

# ╔═╡ 9bb72ff3-9307-4904-a8b7-e3021b64c343
md"""
(예제3) 첫번째 고유값에 대응하는 고유벡터
"""

# ╔═╡ df8db08d-bdda-461b-93cf-e183fc460a41
let 
	A = [1 2 ; 2 1]
	ψ = [-1,1]
	λ = -1 
	A*ψ == λ*ψ
end

# ╔═╡ db0b8cb1-2978-4441-9d5c-5e244f31eb0a
md"""
-  $\psi = \begin{bmatrix} -1 \\ 1 \end{bmatrix}$은 $\lambda=-1$에 대응하는 ${\bf A}$의 고유벡터이다. 
"""

# ╔═╡ 8d7119f8-1ef7-4295-a962-d14b1b624689
md"""
(예제4) 첫번째 고유값에 대응하는 또 다른 고유벡터
"""

# ╔═╡ 2a3b32af-8c6e-495b-9cb6-872d3a64dcd5
let 
	A = [1 2 ; 2 1]
	ψ = [-2,2]
	λ = -1 
	A*ψ == λ*ψ
end

# ╔═╡ 064f317c-a22d-4220-a944-7d2979dc3515
md"""
-  $\psi = \begin{bmatrix} -2 \\ 2 \end{bmatrix}$ 역시 $\lambda=-1$에 대응하는 ${\bf A}$의 고유벡터이다. 
"""

# ╔═╡ d7dba18d-7ee2-41df-b970-d6b92599b22f
md"""
`-` 예제3,4의 관찰: $\psi$가 ${\bf A}$의 고유벡터이라면 $-\psi, \frac{1}{\sqrt{2}}\psi,\dots$ 모두 ${\bf A}$의 고유벡터이다. 그리고 이때 $\psi, -\psi, \frac{1}{\sqrt{2}}\psi$에 대응하는 고유값은 모두 같다. 
"""

# ╔═╡ b85ea93b-70b9-4ba3-894f-8864626fcca1
md"""
(예제5) 두번째 고유값에 대응하는 고유벡터
"""

# ╔═╡ 23b12658-c5af-4b28-a49c-ae6b10815de0
let 
	A = [1 2 ; 2 1]
	ψ = [1,1]
	λ = 3 
	A*ψ == λ*ψ
end

# ╔═╡ 40b14b87-6d3e-47c4-9177-dfeeb116c858
md"""
-  $\psi = \begin{bmatrix} 1 \\ 1 \end{bmatrix}$ 은 $\lambda=3$에 대응하는 ${\bf A}$의 고유벡터이다.
"""

# ╔═╡ 3ede1cfb-87d6-40a7-884c-0da963c24d92
md"""
### 고유값의 존재 
"""

# ╔═╡ 003a2f97-46cf-4465-83f9-0ad4c64b8c59
md"""
`-` **고유값이 없는 정사각행렬은 없다.**

(왜?) 행를 ${\bf A}_{n\times n}$의 고유값이 없다는 의미는 $\det({\bf A}-\lambda {\bf I})=0$를 만족하는 $\lambda$가 없다는 의미이다. 그런데 임의의 $n$차 다항식의 해는 항상 존재한다 (대수학의 기본정리). 
"""

# ╔═╡ fdba1f70-d0c2-4765-91fd-452fdfe487b9
md"""
`-` 그런데 $n$차 다항식의 해가 중복근일 수도 있으므로 ${\bf A}$가 서로 다른 $n$개의 고유값을 가질 필요는 없다.
"""

# ╔═╡ 80393f7b-c7d2-4972-9ed6-6d06439760ee
md"""
(예제1) 2개의 고유값이 겹치는 예제
"""

# ╔═╡ 0757d40b-2402-445e-ad4c-2daaec0b1b75
let 
	A = [1 0; 0 1]
	eigvals(A)
end

# ╔═╡ a9970ddd-eba6-4eb6-9a32-a6bcc154ed68
md"""
(예제2) 2개의 고유값이 겹치는 두번째 예제
"""

# ╔═╡ 719cdcf6-d062-48d3-9371-a6419aa69522
let 
	A = [0 0; 0 0]
	eigvals(A)
end

# ╔═╡ 7e88ea4a-61ff-4fb6-bf6b-14aed6570cce
md"""
### 특성방정식의 근에 대응하는 고유벡터의 존재
"""

# ╔═╡ fa606cc4-7cf2-457f-889a-3a4bdeaaab31
md"""
`-` **특성방정식을 만족하는 근 $\lambda$에 대응하는 고유벡터가 반드시 하나는 존재한다.** 
"""

# ╔═╡ 85302350-4d8f-4a26-8901-b308fd0b7e20
md"""
(왜?) 특성방정식을 만족하는 하나의 근 $\lambda^*$를 fix하자. 고정된 $\lambda^*$에 대응하는 고유벡터가 없다는 의미는 

$$({\bf A}-\lambda^*{\bf I})\psi =0$$

를 만족하는 $\psi$는 오직 $\psi=0$ 뿐이라는 것을 의미한다. 그런데 이는 사실이 아니다. 왜냐하면 $\lambda^*$는 

$$\det({\bf A}-\lambda^*{\bf I})=0$$

를 만족하고 따라서 행렬 ${\bf A}-\lambda^*{\bf I}$는 역행렬이 없는 행렬이 되고 이렇게 되면 ${\bf A}-\lambda^*{\bf I}$의 column 들은 선형독립이 아니게 된다. 따라서 $({\bf A}-\lambda^*{\bf I})\psi=0$을 만족하는 $\psi\neq0$가 있다.
"""

# ╔═╡ 465232b9-61e3-419e-a38e-e26472fa389d
md"""
(보충학습) 벡터 $V_1,V_2,\dots,V_n$이 선형독립이 아니라는 의미는 적당한 $c_1,c_2,\dots,c_n$이 존재하여 

$$c_1 V_1 + c_2V_2 + \dots + c_nV_n=0$$

을 만족한다는 의미이다. (단, 이때 $c_1,\dots c_n$이 모두 0은 아니라고 가정한다.)
"""

# ╔═╡ 3fbd4bb7-d325-4982-ab0a-def127f91434
md"""
- (응용버전) 적당한 매트릭스 ${\bf V}=[V_1~ V_2~ \dots ~ V_n]$에 대하여 ${\bf c} \neq 0$인 벡터 ${\bf c}= (c_1,\dots,c_n)^\top$이 존재하여 $${\bf V}{\bf c}=0$$  만족하면 ${\bf V}$의 column들은 선형독립이 아니라고 볼 수 있다.
"""

# ╔═╡ 2b4c2b61-53b6-4585-9e5e-774a0cb9bc9f
md"""
### 고유벡터의 차원
"""

# ╔═╡ 10a68b69-f319-48ac-acfb-c71964edd2e7
md"""
`-` "하나의 고유값에 반드시 하나의 고유벡터는 존재한다"라는 말은 사실 무한개의 고유벡터가 존재한다는 것을 의미한다. 왜냐하면 $\psi$가 고유벡터이면 $\sqrt{2}\psi$도 고유벡터이기 때문이다. 
"""

# ╔═╡ cac169f2-0fd8-4ca9-9932-4ed751f834bf
md"""
`-` 아래의 예제들을 관찰하자.
"""

# ╔═╡ 23729642-3a1e-4f22-b05b-356fad1bf54a
md"""
(예제1)
"""

# ╔═╡ bc325027-c401-410d-be4e-2e99e5077269
let 
	A = [1 0 
		 0 2]
	eigen(A)
end

# ╔═╡ 73808ba1-63b7-4157-b606-26814a864f2f
md"""
-  $\begin{bmatrix} 1 \\ 0 \end{bmatrix}$ 이 고유벡터 이므로 $\begin{bmatrix} 0.5 \\ 0 \end{bmatrix}$, $\begin{bmatrix} -3.14 \\ 0 \end{bmatrix}$ 등도 고유벡터이다. 
"""

# ╔═╡ 2243a0ee-3329-4b06-883c-6dd8bfd44551
md"""
-  $\begin{bmatrix} 0 \\ 1 \end{bmatrix}$ 이 고유벡터 이므로 $\begin{bmatrix} 0 \\ 5.3 \end{bmatrix}$, $\begin{bmatrix} 0 \\ -\sqrt{3} \end{bmatrix}$ 등도 고유벡터이다. 
"""

# ╔═╡ 14c834bd-61ef-45d8-8aa1-cfee057f1e3d
md"""
(예제2)
"""

# ╔═╡ cbdabff5-c038-4f7a-ac35-545aa7127bd3
let 
	B = [1 1 
		 0 1]
	eigen(B)
end

# ╔═╡ 212b9778-dee0-45a2-a4a8-b75f4fa7a0bb
md"""
-  $\begin{bmatrix} 1 \\ 0 \end{bmatrix}$ 이 고유벡터 이므로 $\begin{bmatrix} -1 \\ 0 \end{bmatrix}$, $\begin{bmatrix} -3.14 \\ 0 \end{bmatrix}$ 등도 고유벡터이다. 
"""

# ╔═╡ 8657dff6-92e6-45df-853d-d13ee2959138
md"""
`-` 예제1,2 모두 무한개의 고유벡터를 가지는 것은 맞지만 차원이 다르다. 예제2의 경우 $\begin{bmatrix} 1 \\ 0\end{bmatrix}$을 basis로 한 벡터들의 조합을 만들 수 있지만 예제1의 경우 $\begin{bmatrix} 1 \\ 0\end{bmatrix}$을 basis로 한 벡터들의 조합을 만들 수 있고 추가적으로 $\begin{bmatrix} 0 \\ 1\end{bmatrix}$를 basis로 한 벡터들의 조합도 만들 수 있기 때문이다. 
"""

# ╔═╡ fa6d45e7-3921-4b86-a757-1c1eddda5438
md"""
`-` 예제1이 가질 수 있는 고유벡터 조합들이 예제2보다 더 풍부하다. 이 느낌을 좀 더 수학적으로 표현할 수 있을까? $\to$ 예제1의 고유벡터들로 조합할 수 있는 (=고유벡터들의 조합으로 확장할 수 있는, 고유벡터로 span하는) 점들은 2차원 평면을 만들지만 예제2의 고유벡터로 조합할 수 있는 (=고유벡터들의 조합으로 확장할 수 있는, 고유벡터로 span하는) 점들은 1차원 직선을 만든다. 
"""

# ╔═╡ f32f1283-1689-4372-b89b-697604efe6fc
md"""
(예제3) 예제1,2의 고유벡터로 조합할 수 있는 (=고유벡터들의 조합으로 확장할 수 있는, 고유벡터로 span하는) 점들을 시각화 하라. 예제1의 매트릭스를 임의로 바꿔보면서 관찰해보라. 
"""

# ╔═╡ 71dbec22-42b4-4cbf-9e86-e8952ec127c4
md"c1 $(@bind c1 Slider(-5:0.1:5,show_value=true))"

# ╔═╡ 95674f32-720c-48bc-8d88-f8ed33b270ae
md"c2 $(@bind c2 Slider(-5:0.1:5,show_value=true))"

# ╔═╡ 0a987b69-88dc-4ba2-af1c-c39cb5d85019
let
	A = [1 5 
		 2 2]
	B = [1 1 
		 0 1]
	f = Ψ -> [Ψ[:,i] for i in 1:2]
	v1, v2= A |> eigvecs |> f
	u1, u2= B |> eigvecs |> f
	Bdx, Bdy = [-5,5,-5,5]*v1' + [5,5,-5,-5]*v2' |> f
	Ax,Ay = c1*v1 + c2*v2
	Bx,By = c1*u1 + c2*u2
	scatter(xlim=(-11,11),ylim=(-11,11))
	scatter!(Bdx,Bdy,markershape=:cross,alpha=0.2,color=:"red")
	scatter!([Ax],[Ay],markershape=:cross,color=:"red")
	scatter!([Bx],[By],color=:"blue")
end 

# ╔═╡ de97e015-f2d8-453e-9f71-e748eab4159e
md"""
## 숙제 

정사각행렬 ${\bf A}_{2\times 2}$가 0행렬인 경우에 ${\bf A}$의 고유벡터들이 span하는 공간이 몇차원인가? 
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Plots = "~1.29.0"
PlutoUI = "~0.7.39"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9489214b993cd42d17f44c36e359bf6a7c919abf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "7297381ccb5df764549818d9a7d57e45f1057d30"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.18.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "a985dc37e357a3b22b260a5def99f3530fb415d3"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.2"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "b153278a25dd42c65abbf4e62344f9d22e59191b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.43.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "cc1a8e22627f33c789ab60b36a9132ac050bbf75"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.12"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "51d2dfe8e590fbd74e7a842cf6d13d8a2f45dc01"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.6+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "b316fd18f5bc025fedcb708332aecb3e13b9b453"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.3"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "1e5490a51b4e9d07e8b04836f6008f46b48aaa87"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.3+0"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "336cc738f03e069ef2cac55a104eb823455dca75"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.4"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "46a39b9c58749eefb5f2dc1178cb8fab5332b1ab"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.15"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "d457f881ea56bbfa18222642de51e0abf67b9027"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.29.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "dc1e451e15d90347a7decc4221842a022b011714"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.5.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "cdbd3b1338c72ce29d9584fdbe9e9b70eeb5adca"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.1.3"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "bc40f042cfcc56230f781d92db71f0e21496dffd"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.5"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "cd56bf18ed715e8b09f06ef8c6b781e6cdc49911"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c82aaa13b44ea00134f8c9c89819477bd3986ecd"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.3.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "e75d82493681dfd884a357952bbd7ab0608e1dc3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.7"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─e08c8300-db18-11ec-3f60-19c9638fa2d2
# ╟─90db85a9-9864-4e63-98ca-7edb77f1932f
# ╟─c84fb7f0-e572-46ab-b0a3-7612c1870815
# ╠═15442a71-00fe-48ed-ae51-def2f2d8876e
# ╠═8308a4c4-46b0-4407-aba5-d8b174c2f152
# ╠═8b5ce365-94e4-4875-ac48-d7d6f4b7409c
# ╟─9d9020be-5c18-48f6-bf87-d3623f568673
# ╟─05c1b786-159f-4ebd-9748-65c5194133f6
# ╟─65396179-2be2-4ba2-84a6-38db6a4335c5
# ╟─318c6975-764e-47e5-abe7-0ba7b9883d20
# ╠═a411e2b7-b659-44d1-8146-16d15f59e0e4
# ╟─867a6ec4-fb9c-4001-910f-9080ecfcddbc
# ╟─67b61ca4-a06b-4f1d-8cab-0d4ba12dfae6
# ╠═11d55eb6-67af-4146-abe0-c34dba602e71
# ╟─6f9a1a8e-717b-4953-804a-e8647743fbd4
# ╟─598d1691-16cf-4803-a5af-d87abec1a02e
# ╟─fa178a75-8119-4d2c-9f63-94420591e7a5
# ╟─c1e1aded-6005-4b01-bab1-8f75a3811dc4
# ╟─752ebf31-3ec1-4033-a3d1-0bddd8ac3232
# ╟─3a85c442-10cd-46ac-943f-6912018b008c
# ╟─330e2f3a-4d66-417f-aaf5-b92ba0a79c55
# ╟─713775f6-ccce-4b7a-bda3-e828aebcd6a1
# ╟─006740e0-f12c-4312-b181-30d5e6f82ec3
# ╟─cba26403-d491-4eae-be95-7a11984c2b0e
# ╟─fe5876ba-30f9-470d-bd42-4c188fa580f3
# ╟─bca97f96-16d9-4c2d-8d6b-6bd5c7400706
# ╠═ee444ee9-9a95-4b8f-999d-60d1c16df56f
# ╟─4a687148-264d-46b4-a1e4-7d3e94288e57
# ╟─17bce861-3ef6-4050-bc79-67da0ea689b8
# ╠═6349ab90-43fd-4a19-87d2-aa2483d3faea
# ╟─d7d62756-2dab-4a0f-817a-e1acd83a3b74
# ╟─9bb72ff3-9307-4904-a8b7-e3021b64c343
# ╠═df8db08d-bdda-461b-93cf-e183fc460a41
# ╟─db0b8cb1-2978-4441-9d5c-5e244f31eb0a
# ╟─8d7119f8-1ef7-4295-a962-d14b1b624689
# ╠═2a3b32af-8c6e-495b-9cb6-872d3a64dcd5
# ╟─064f317c-a22d-4220-a944-7d2979dc3515
# ╟─d7dba18d-7ee2-41df-b970-d6b92599b22f
# ╟─b85ea93b-70b9-4ba3-894f-8864626fcca1
# ╠═23b12658-c5af-4b28-a49c-ae6b10815de0
# ╟─40b14b87-6d3e-47c4-9177-dfeeb116c858
# ╟─3ede1cfb-87d6-40a7-884c-0da963c24d92
# ╟─003a2f97-46cf-4465-83f9-0ad4c64b8c59
# ╟─fdba1f70-d0c2-4765-91fd-452fdfe487b9
# ╟─80393f7b-c7d2-4972-9ed6-6d06439760ee
# ╠═0757d40b-2402-445e-ad4c-2daaec0b1b75
# ╟─a9970ddd-eba6-4eb6-9a32-a6bcc154ed68
# ╠═719cdcf6-d062-48d3-9371-a6419aa69522
# ╟─7e88ea4a-61ff-4fb6-bf6b-14aed6570cce
# ╟─fa606cc4-7cf2-457f-889a-3a4bdeaaab31
# ╟─85302350-4d8f-4a26-8901-b308fd0b7e20
# ╟─465232b9-61e3-419e-a38e-e26472fa389d
# ╟─3fbd4bb7-d325-4982-ab0a-def127f91434
# ╟─2b4c2b61-53b6-4585-9e5e-774a0cb9bc9f
# ╟─10a68b69-f319-48ac-acfb-c71964edd2e7
# ╟─cac169f2-0fd8-4ca9-9932-4ed751f834bf
# ╟─23729642-3a1e-4f22-b05b-356fad1bf54a
# ╠═bc325027-c401-410d-be4e-2e99e5077269
# ╟─73808ba1-63b7-4157-b606-26814a864f2f
# ╟─2243a0ee-3329-4b06-883c-6dd8bfd44551
# ╟─14c834bd-61ef-45d8-8aa1-cfee057f1e3d
# ╠═cbdabff5-c038-4f7a-ac35-545aa7127bd3
# ╟─212b9778-dee0-45a2-a4a8-b75f4fa7a0bb
# ╟─8657dff6-92e6-45df-853d-d13ee2959138
# ╟─fa6d45e7-3921-4b86-a757-1c1eddda5438
# ╟─f32f1283-1689-4372-b89b-697604efe6fc
# ╠═71dbec22-42b4-4cbf-9e86-e8952ec127c4
# ╠═95674f32-720c-48bc-8d88-f8ed33b270ae
# ╠═0a987b69-88dc-4ba2-af1c-c39cb5d85019
# ╟─de97e015-f2d8-453e-9f71-e748eab4159e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
