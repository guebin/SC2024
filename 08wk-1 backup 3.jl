### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 3c5aa44e-d594-11ec-21d6-25d65911031d
using PlutoUI,Plots, Distributions, LinearAlgebra,RDatasets

# ╔═╡ 1f42323f-5074-41f8-bd6a-1806c123e149
md"""
# 08wk-1: 주성분분석 
"""

# ╔═╡ 67efeeff-8a5c-476f-ab22-f3ff4208439c
md"""
## 1. 강의영상
"""

# ╔═╡ 32bd60e5-2872-430a-8053-421b6e454f83
html"""
<div style="display: flex; justify-content: center;">
<div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
<iframe src=
"
https://youtube.com/embed/playlist?list=PLQqh36zP38-zRx2DBl6k4gwIqRNdjrRYf&si=tJxKKXnNxHE9bvJ_
"
width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""
	

# ╔═╡ e36c494d-0e97-4010-a57d-fcd60ab93900
md"""
## 2. Imports
"""

# ╔═╡ 107c213d-b88e-4138-9a1e-0fba804f6e09
PlutoUI.TableOfContents()

# ╔═╡ 9af6ddd2-d726-4ed0-baa1-0c692d1f96d9
Plots.plotly()

# ╔═╡ ece414e2-e7dd-40f5-b401-5f5039d9b804
md"""
## 3. SVD의 응용: 주성분분석
"""

# ╔═╡ b3ab776c-9dbc-4e17-9f28-0c8f26bba169
md"""
-- 주성분분석(PCA)은 차원축소에 이용되거나 다중공선성의 해결에 이용되는 기법이다. 
"""

# ╔═╡ 174ef9ba-a3c0-4040-8d31-1cd6962b7b8f
md"""
-- 우선은 차원축소의 관점에서만 PCA를 서술해보자. 
"""

# ╔═╡ d56fac91-c4c1-4529-9a75-8f4e007cf258
md"""
### A. PCA의 목적
"""

# ╔═╡ 45ca0b83-4ca4-4137-9cbf-a03ca32d4897
md"""
-- 데이터: 아래와 같은 $n\times p$ 매트릭스 (혹은 데이터프레임)이 있다고 하자. 

${\bf X}_{n\times p}$

- 통계학과에서 많이 쓰는 그 design matrix라고 생각하면 된다. 
- 여기에서 $n>p$ 를 가정한다. 

"""

# ╔═╡ 762c7cf9-73ab-4fce-9f74-177140560eec
md"""
-- 소망: `(1)` ${\bf X}$가 가지고 있는 정보는 거의 유지하면서 `(2)` 저장비용은 좀 더 저렴한 매트릭스 ${\bf Z}$를 만들고 싶다. (얌체아냐??)
- 정보를 거의 유지한다는 것이 무슨 말? 
- 저장비용이 저렴하다는 것은 무슨 말?

"""

# ╔═╡ 989e6e04-1d01-41be-9d77-973089e0264b
md"""
-- 소망을 좀 더 구체적으로 수식화하면 아래와 같다. 

`(1)` ${\bf Z}$의 차원이 $(n,q)$ 이며, (단 여기에서 $q<p$). 

`(2)` ${\bf Z}$에서 적당한 변환 ${\bf B}$를 하면 언제든지 ${\bf X}$와 비슷한 매트릭스로 복원(reconstruction) 할 수 있다. 즉

${\bf Z}{\bf B} \approx {\bf X}$

을 만족하는 적당한 ${\bf B}$가 존재한다. (이렇게 되면 변환 ${\bf B}$ 의 치원은 $(q,p)$ 가 되어야 한다.)

"""

# ╔═╡ 9c968667-0a77-42ce-b381-8aa3cf5e31e2
md"""
### B. SVD $\to$ ${\bf Z}$를 찾아보자.
"""

# ╔═╡ 0c36fb86-88e4-4f92-8b34-94944150f273
md"""
-- 우선 임의의 매트릭스 ${\bf X}_{n\times p}$의 SVD를 이용하면 적당히 아래를 만족하는 매트릭스들을 찾을 수 있다. 

- ``{\bf X}={\bf U}{\bf D}{\bf V}^\top= {\bf U}_1 {\bf D}_1 {\bf V}_1^\top+ {\bf U}_2 {\bf D}_2 {\bf V}_2^\top``
- ``{\bf X} \approx {\bf U}_1 {\bf D}_1 {\bf V}_1^\top`` *<-- 이렇게 잡을수 있는게 포인트!!*

여기에서 각 매트릭스의 차원은 

- ``{\bf U}_1``: $(n,q)$
- ``{\bf U}_2``: $(n,p-q)$
- ``{\bf D}_1``: $(q,q)$
- ``{\bf D}_2``: $(p-q,p-q)$
- ``{\bf V}_1``: $(p,q)$, ``{\bf V}_1^\top = (q,p)``
- ``{\bf V}_2``: $(p,p-q)$, ``{\bf V}_2^\top = (p-q,p)``
"""

# ╔═╡ ffd6d9f2-e4a5-4f06-a75b-242cb4ee71d1
md"""
-- 여기에서 ${\bf U}_1{\bf D}_1$를 유심히 관찰해보자. 아래와 같은 특징이 있다. 

`(1)` ${\bf U}_1{\bf D}_1$의 차원이 $(n,q)$ 이다. 

`(2)` ${\bf U}_1{\bf D}_1$에서 ${\bf V}_1^\top$를 곱하면 ${\bf X}$를 비슷하게 만들 수 있다.  ($\because {\bf U}_1{\bf D}_1{\bf V}_1 \approx {\bf X}$)

따라서 

${\bf Z}:={\bf U}_1{\bf D}_1, \quad {\bf B}:={\bf V}^\top$

와 같이 정의하면 PCA의 소망에서 언급한 `(1)`,`(2)`의 조건을 만족하게 된다. 
"""

# ╔═╡ 0ae3053c-51de-4fd0-809d-0df5f539de61
md"""
!!! info "주성분 분석"
	임의의 $(n,p)$ matrix ${\bf X}$ 에 대한 SVD가 아래와 같다고 하자. $(n \geq p)$.
	
	${\bf X} = {\bf U}{\bf D}{\bf V}^\top$

	이때 ${\bf Z}={\bf U}{\bf D}={\bf X}{\bf V}$ 의 each-column 을 주성분 (principle componet) 라고 한다. 구체적으로 
	- PC1 = first column of ${\bf Z}$. 
	- PC2 = second column of ${\bf Z}$. 
	- ``\dots``
	- PCp = $p$-th column of ${\bf Z}$. 
	와 같이 부른다. 또한 ${\bf Z}$ 의 각 원소를 principal component 라고 부른다. 따라서 ${\bf Z}_{n \times p}$ 에는 $np$ 개의 principal component score 가 있다. 그리고 ${\bf V}$ 는 rotation matrix 혹은 loading matrix 라고 부른다. 마지막으로 원래의 자료 ${\bf X}$를 그대로 분석하것이 아니라, ${\bf X}$를 ${\bf Z}$로 바꾼뒤에 분석하는 일련의 기법 (즉 ${\bf X}$의 주성분을 분석하는 기법)을 통칭하여 주성분분석이라고 한다. 
"""

# ╔═╡ dd3262c1-d3fb-4feb-ac5b-187833541b19
md"""
### C. IRIS data
"""

# ╔═╡ cee6edd8-a040-4c7c-a417-854904fb413b
md"""
-- iris 데이터 로드
"""

# ╔═╡ ba4fdd2e-f83b-4b85-9ec2-85119c97ab15
iris = dataset("datasets","iris")

# ╔═╡ 0017fd0e-2ba5-4c30-a2f1-b0f6edf0c3c2
md"""
-- 각 열에 접근
"""

# ╔═╡ c9d2ece5-6594-454c-9e7a-20a430a08751
iris.SepalLength # 방법1

# ╔═╡ f5acceb2-3523-4f66-ae54-4ea8fb01ee12
iris[:,1] # 방법2

# ╔═╡ 48178775-ed6f-45e8-bfcf-dbc1a93fa98d
md"""
-- 데이터의 정리
"""

# ╔═╡ 3620eaa8-5ed5-45c7-bac1-aaa2de363de4
X,y = Array(iris[:,1:4]), iris[:,5]

# ╔═╡ 238de6fe-735d-4ada-9543-3b78ff3f4b2a
md"""
-- SVD를 이용한 주성분 분석 수행
"""

# ╔═╡ 27548de8-5e9a-407f-b3e4-97b596046dfb
let 
	U,d,V = svd(X) # (150,4)
	U1,U2 = U[:,1:2],U[:,3:4]
	D1,D2 = Diagonal(d[1:2]), Diagonal(d[3:4])
	V1,V2 = V[:,1:2],V[:,3:4]
	Z = U1*D1 # (150,2)
	#[X Z*V1']
end 

# ╔═╡ 92108017-9d7b-4434-be5f-1dd3c3811131
md"""
-- 주성분분석을 이용한 시각화
"""

# ╔═╡ 759ffa71-0217-4baf-8e12-084c7bb004c5
let
	U,d,V = svd(X) # (150,4)
	U1,U2,U3,U4 = eachcol(U)
	d1,d2,d3,d4 = d
	V1,V2,V3,V4 = eachcol(V)
	#Z1,Z2 = U1*d1, U2*d2
	PC1,PC2 = U1*d1, U2*d2
	@show unique(y)
	scatter(PC1[y .== "setosa"],PC2[y .== "setosa"],label="setosa")
	scatter!(PC1[y .== "versicolor"],PC2[y .== "versicolor"],label="versicolor")
	scatter!(PC1[y .== "virginica"],PC2[y .== "virginica"],label="virginica")
end 

# ╔═╡ 00367152-c595-4a2e-9666-01d1a136b4f8
md"""
-- Perfect Reconstruction
"""

# ╔═╡ 8f931cbf-f2fc-41e2-992a-3061756b558c
let
	U,d,V = svd(X)
	Z = U*Diagonal(d)
	Z*V'
end 

# ╔═╡ 502c4e40-520b-4051-92d8-243e97ed7982
md"""
### D. PCA에 대한 자잘이들.
"""

# ╔═╡ 839ee40e-8e0a-409d-ac02-fad2233b7454
md"""
!!! info "자잘이1: 주성분분석시 보통 scaling을 선행한다"
	주성분분석을 할때 ${\bf X}$ 자체를 특이값분해하기 보다 ${\bf X}$를 열별로 표준화한것을 고려하는 경우가 많다. (이유는 ${\bf X}$가 가지고 있는 변수사이의 관계를 포착하기 위해서!!)
"""

# ╔═╡ 14ee9809-e4dd-41bc-ae62-669606e954d0
let 
	X1,X2,X3,X4 = eachcol(X)
	X1 = (X1 .- mean(X1))/std(X1)
	X2 = (X2 .- mean(X2))/std(X2)
	X3 = (X3 .- mean(X3))/std(X3)
	X4 = (X4 .- mean(X4))/std(X4)
	[X1 X2 X3 X4] # 보통 이걸 X라고 생각하고 분석함
end 

# ╔═╡ 3e64482f-7fc6-4544-ab28-189dc7293830
md"""
!!! info "자잘이2: PC score를 쓰는 경우"
	PC score는 보통 observation 별로 분석하고 싶을 때 사용한다. 따라서 **"PC scores of $i$-th observation**" 와 같은 문맥으로 자주 쓴다. 
"""

# ╔═╡ fe973b50-53a8-49c6-bd67-ebd947ca1cc9
let 
	U,d,V = svd(X)
	Z = U * Diagonal(d)
	Z[1,:] # PC scores of first obs 
end 

# ╔═╡ cbe9a2f9-2874-4635-bac9-3b51d1a50a78
md"""
!!! info "자잘이3: ?? 개의 주성분이라는 표현"
	``{\bf X}_{n \times p}`` 인 행렬에서 주성분의 숫자는 최대 $p$개 까지 있을 수 있다. 만약에 **"${\bf X}$의 주성분을 고려하면"** 이라는 용어를 쓴다면 아래와 같은 ${\bf Z}$를 고려한다는 의미이다. 
	
	```julia
	U,d,V = svd(X) # X: (n,p)
	Z = U*Diagonal(d) # Z: (n,p)
	```

	만약에 **"${\bf X}$의 (처음) 3개 주성분만 고려하면"** 이라는 용어를 쓴다면 아래와 같은 ${\bf Z}$를 고려한다는 의미이다. 

	```julia
	U,d,V = svd(X) # X: (n,p)
	Z = U[:,1:3]*Diagonal(d[1:3]) # Z: (n,3)
	```
"""

# ╔═╡ 2e2a7afe-b2ff-4683-9ab0-8474e946e8be
md"""
## 4. ${\bf Z}$를 계산하는 여러 방법
"""

# ╔═╡ b7ead0c1-079d-48ee-bf4c-aa699bb50161
md"""
### A. ${\bf Z} = {\bf UD}$
"""

# ╔═╡ 4ef620f5-de71-492c-bf26-e979d9c1ff54
let 
	U,d,V = svd(X)
	Z = U*Diagonal(d)
end 

# ╔═╡ 27238d09-b7d9-46ad-b33e-937f1ffe449e
md"""
### B. ${\bf Z} = {\bf X} {\bf V}$
"""

# ╔═╡ a903457f-43a5-49b6-9d28-7eaff4c67ced
let 
	U,d,V = svd(X)
	Z = X*V 
end 

# ╔═╡ 62a9ea09-ae87-40e7-893b-1a5449100c13
md"""
### C. ${\bf X}^\top{\bf X}$의 SVD를 이용
"""

# ╔═╡ 56dc1c61-22d3-4f10-b5df-5869b2cff319
let 
	U,d,V = svd(X'X)
	@show V ≈ U
	X*V
	#X*U
end 

# ╔═╡ d27ab23b-efd4-461d-9266-43c280c64dd6
md"""
### D. ${\bf Z} = {\bf X \Psi}$
"""

# ╔═╡ d7ef4f20-a885-4242-bb60-b2fad3d6b511
let 
	λ,Ψ = eigen(X'X, sortby = -)
	X*Ψ
end 

# ╔═╡ 0239f3aa-2db4-496c-9920-b9d1704afe8a
let 
	U,d,V = svd(X'X)
	λ,Ψ = eigen(X'X, sortby = -)
	@show U ≈ V ≈ Ψ
	@show d ≈ λ
end 

# ╔═╡ fcb48c59-e26c-4c17-b087-beab10b3b66e
md"""
## 5. "고유값분해 = SVD" 인 경우
"""

# ╔═╡ 6fed3d84-bb66-4951-9dae-ddc221e14669
md"""
!!! info "'고유값분해 = 특이값분해' 인 경우"
	정사각행렬 $\bf X$ 가 어떠한 조건을 만족하면 특이값분해의 ${\bf U}$와 ${\bf V}$가 같아지는 경우가 있다. 즉 ${\bf X}$가 아래와 같이 분해되는 경우가 있다. 

	${\bf X}={\bf V}{\bf D}{\bf V}^\top={\bf U}{\bf D}{\bf U}^\top$

	이 상황을 코드로 표현하면 아래와 같다. 

	```julia
	U,d,V = svd(X)
	U ≈ V # <- true
	```

	또한 이때 ${\bf V}$는 (매우 놀랍게도) 아래를 실행하여 얻어지는 고유벡터 행렬 ${\bf \Psi}$와 일치하게 된다.

	```julia
	λ, Ψ = eigen(X, sortby = -)
	```

	따라서 이 내용을 종합하면 ${\bf X}$가 특정한 조건을 만족할 때, 아래가 성립하는 경우가 있다고 기억하면 된다. 
	
	```julia
	λ, Ψ = eigen(X, sortby = -)
	U,d,V = svd(X)
	U ≈ V ≈ Ψ # true
	d ≈ λ # true
	```

	수식으로 쓰면 아래와 같이 쓸 수 있다. 

	${\bf X}={\bf U}{\bf D}{\bf V}^\top={\bf \Psi}{\bf \Lambda}{\bf \Psi}^\top$

	그렇다면 이러한 경우를 만드는 "특수한 ${\bf X}$" 는 무엇일까?
"""

# ╔═╡ 4c934d1d-7962-4076-a8ce-a1c7b2788652
md"""
### A. (실대칭행렬 + $\alpha$) 의 조건에서 가능
"""

# ╔═╡ 1c8e324f-c9b0-40cf-862c-6e9f726d6e0d
md"""
-- 예제1: 실대칭행렬이 아니어서 불가능한 경우
"""

# ╔═╡ dfc79091-f1c3-4ca4-812b-bd3b95bff483
# 대칭행렬이 아니야..
let
	X = [0 1 0 
		 0 0 1 
	     1 0 0]
	U,d,V = svd(X)
	@show U ≈ V 
end 

# ╔═╡ 2789d231-f61d-469d-9ba8-d65afd99da0e
md"""
-- 예제2: 실대칭행렬이 아니어서 불가능한 경우2
"""

# ╔═╡ 37750ce5-0c0c-488a-9e35-4274b4860d67
# 대칭행렬이지만 원소가 실수가 아니야..
let
	X = [2im 1
		 1   0]
	U,d,V = svd(X)
	@show U ≈ V
end 

# ╔═╡ 14b7545c-9c63-449f-9b01-b7df40914694
md"""
-- 예제3: 실대칭행렬이지만 불가능한 경우
"""

# ╔═╡ 21038d5f-9017-4860-ac1f-c4ca7c8f8037
# 실대칭행렬이지만 불가능함
let 
	X = [2  0
		 0 -1]
	U,d,V = svd(X)
	@show U ≈ V
end 

# ╔═╡ d4bff255-c5c8-43c1-9493-2184deaab5ad
md"""
-- 예제4: (실대칭행렬+$\alpha$)의 조건에서 가능한 경우
"""

# ╔═╡ 28a12ec6-ac27-429c-b183-de313b359100
# 실대칭행렬 + alpha 의 조건에서 가능
let 
	X = [1 2
	     2 4]
	U,d,V = svd(X)
	@show U ≈ V
	λ, Ψ = eigen(X, sortby = -)
	U = V = Ψ # 이렇게 U,V 를 "선택"하자..
	@show X ≈ U * Diagonal(d) * V' # -- SVD의 정의만족
	@show X ≈ Ψ * Diagonal(λ) * Ψ' # -- Eigen Decomposition의 정의만족
end 

# ╔═╡ b4e47ad5-656d-406c-9268-d95c5e59028f
md"""
### B. 무조건 가능한 경우
"""

# ╔═╡ 7f6defd4-670e-4282-ab67-81ff268e7a6b
md"""
*통계학과에서 자주 사용하는 design matrix ${\bf X}_{n\times p}$, $n \geq p$를 생각하자.*
"""

# ╔═╡ 2967744e-7512-4b87-9811-9c35e4813590
md"""
-- 예제1: 공분산행렬
"""

# ╔═╡ f57616cd-26c8-4e68-9c4e-91236b10718e
let
	U,d,V = svd(cov(X))
	@show U ≈ V
	λ,Ψ = eigen(cov(X),sortby= - )
	@show U ≈ V ≈ Ψ
	U = V = Ψ # 이렇게 U,V 를 "선택"하자..
	@show cov(X) ≈ U * Diagonal(d) * V' # -- SVD의 정의만족
	@show cov(X) ≈ Ψ * Diagonal(λ) * Ψ' # -- Eigen Decomposition의 정의만족
end 

# ╔═╡ 2e2bfeb9-66f5-4d92-8a07-3256b5dcb742
md"""
-- 예제2: 상관계수 행렬
"""

# ╔═╡ ab13954d-5c71-4fde-9e33-11016ca8368d
let
	U,d,V = svd(cor(X))
	@show U ≈ V
	λ,Ψ = eigen(cor(X),sortby= - )
	@show U ≈ V ≈ Ψ
	U = V = Ψ # 이렇게 U,V 를 "선택"하자..
	@show cor(X) ≈ U * Diagonal(d) * V' # -- SVD의 정의만족
	@show cor(X) ≈ Ψ * Diagonal(λ) * Ψ' # -- Eigen Decomposition의 정의만족
end 

# ╔═╡ b9a1c683-5b09-4217-9294-c79549b69e6e
md"""
-- 예제3: 그냥 ${\bf X}^\top {\bf X}$의 꼴..
"""

# ╔═╡ 2f17a7f4-68b0-49ee-b744-352720788178
let
	U,d,V = svd(X'X)
	λ,Ψ = eigen(X'X,sortby= - )
	@show U ≈ V ≈ Ψ # 선택하지 않아도 알아서 잘되네 (운좋음)
	@show λ ≈ d
end 

# ╔═╡ 08698ea9-7ca0-41f6-a5be-aaefd0255828
md"""
-- 예제4: hat-matrix
"""

# ╔═╡ ab34e79e-4c3c-49e3-8186-8576c6e0acff
let 
	H = X*inv(X'X)X'
	U,d,V = svd(H)
	@show U ≈ V # 된다면서요..?
	U1,U2 = U[:,1:4], U[:,5:150]
	d1,d2 = d[1:4], d[5:150]
	V1,V2 = V[:,1:4], V[:,5:150]
	@show U1 ≈ V1 # 여기까지는 잘 됨. 
	@show H ≈ U1*Diagonal(d1)*V1' # U2*Diagonal(d2)*V2' 는 버려도 무방
	@show H ≈ [U1 U2] * Diagonal(d) * [V1 V2]' # 줄리아에서 구해준 V -- SVD 정의만족 
	@show H ≈ [U1 U2] * Diagonal(d) * [V1 U2]' # 내가 내맘대로 바꿔버린 V -- SVD 정의만족
end 

# ╔═╡ 3319fff4-3538-4fca-9210-2a0904c62dff
md"""
- SVD의 분해결과는 유일하지 않음
"""

# ╔═╡ Cell order:
# ╟─1f42323f-5074-41f8-bd6a-1806c123e149
# ╟─67efeeff-8a5c-476f-ab22-f3ff4208439c
# ╟─32bd60e5-2872-430a-8053-421b6e454f83
# ╟─e36c494d-0e97-4010-a57d-fcd60ab93900
# ╠═3c5aa44e-d594-11ec-21d6-25d65911031d
# ╠═107c213d-b88e-4138-9a1e-0fba804f6e09
# ╠═9af6ddd2-d726-4ed0-baa1-0c692d1f96d9
# ╟─ece414e2-e7dd-40f5-b401-5f5039d9b804
# ╟─b3ab776c-9dbc-4e17-9f28-0c8f26bba169
# ╟─174ef9ba-a3c0-4040-8d31-1cd6962b7b8f
# ╟─d56fac91-c4c1-4529-9a75-8f4e007cf258
# ╟─45ca0b83-4ca4-4137-9cbf-a03ca32d4897
# ╟─762c7cf9-73ab-4fce-9f74-177140560eec
# ╟─989e6e04-1d01-41be-9d77-973089e0264b
# ╟─9c968667-0a77-42ce-b381-8aa3cf5e31e2
# ╟─0c36fb86-88e4-4f92-8b34-94944150f273
# ╟─ffd6d9f2-e4a5-4f06-a75b-242cb4ee71d1
# ╟─0ae3053c-51de-4fd0-809d-0df5f539de61
# ╟─dd3262c1-d3fb-4feb-ac5b-187833541b19
# ╟─cee6edd8-a040-4c7c-a417-854904fb413b
# ╠═ba4fdd2e-f83b-4b85-9ec2-85119c97ab15
# ╟─0017fd0e-2ba5-4c30-a2f1-b0f6edf0c3c2
# ╠═c9d2ece5-6594-454c-9e7a-20a430a08751
# ╠═f5acceb2-3523-4f66-ae54-4ea8fb01ee12
# ╟─48178775-ed6f-45e8-bfcf-dbc1a93fa98d
# ╠═3620eaa8-5ed5-45c7-bac1-aaa2de363de4
# ╟─238de6fe-735d-4ada-9543-3b78ff3f4b2a
# ╠═27548de8-5e9a-407f-b3e4-97b596046dfb
# ╟─92108017-9d7b-4434-be5f-1dd3c3811131
# ╠═759ffa71-0217-4baf-8e12-084c7bb004c5
# ╟─00367152-c595-4a2e-9666-01d1a136b4f8
# ╠═8f931cbf-f2fc-41e2-992a-3061756b558c
# ╟─502c4e40-520b-4051-92d8-243e97ed7982
# ╟─839ee40e-8e0a-409d-ac02-fad2233b7454
# ╠═14ee9809-e4dd-41bc-ae62-669606e954d0
# ╟─3e64482f-7fc6-4544-ab28-189dc7293830
# ╠═fe973b50-53a8-49c6-bd67-ebd947ca1cc9
# ╟─cbe9a2f9-2874-4635-bac9-3b51d1a50a78
# ╟─2e2a7afe-b2ff-4683-9ab0-8474e946e8be
# ╟─b7ead0c1-079d-48ee-bf4c-aa699bb50161
# ╠═4ef620f5-de71-492c-bf26-e979d9c1ff54
# ╟─27238d09-b7d9-46ad-b33e-937f1ffe449e
# ╠═a903457f-43a5-49b6-9d28-7eaff4c67ced
# ╟─62a9ea09-ae87-40e7-893b-1a5449100c13
# ╠═56dc1c61-22d3-4f10-b5df-5869b2cff319
# ╟─d27ab23b-efd4-461d-9266-43c280c64dd6
# ╠═d7ef4f20-a885-4242-bb60-b2fad3d6b511
# ╠═0239f3aa-2db4-496c-9920-b9d1704afe8a
# ╟─fcb48c59-e26c-4c17-b087-beab10b3b66e
# ╟─6fed3d84-bb66-4951-9dae-ddc221e14669
# ╟─4c934d1d-7962-4076-a8ce-a1c7b2788652
# ╟─1c8e324f-c9b0-40cf-862c-6e9f726d6e0d
# ╠═dfc79091-f1c3-4ca4-812b-bd3b95bff483
# ╟─2789d231-f61d-469d-9ba8-d65afd99da0e
# ╠═37750ce5-0c0c-488a-9e35-4274b4860d67
# ╟─14b7545c-9c63-449f-9b01-b7df40914694
# ╠═21038d5f-9017-4860-ac1f-c4ca7c8f8037
# ╟─d4bff255-c5c8-43c1-9493-2184deaab5ad
# ╠═28a12ec6-ac27-429c-b183-de313b359100
# ╟─b4e47ad5-656d-406c-9268-d95c5e59028f
# ╟─7f6defd4-670e-4282-ab67-81ff268e7a6b
# ╟─2967744e-7512-4b87-9811-9c35e4813590
# ╠═f57616cd-26c8-4e68-9c4e-91236b10718e
# ╟─2e2bfeb9-66f5-4d92-8a07-3256b5dcb742
# ╠═ab13954d-5c71-4fde-9e33-11016ca8368d
# ╟─b9a1c683-5b09-4217-9294-c79549b69e6e
# ╠═2f17a7f4-68b0-49ee-b744-352720788178
# ╟─08698ea9-7ca0-41f6-a5be-aaefd0255828
# ╠═ab34e79e-4c3c-49e3-8186-8576c6e0acff
# ╟─3319fff4-3538-4fca-9210-2a0904c62dff
