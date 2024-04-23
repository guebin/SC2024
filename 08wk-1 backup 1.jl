### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 3c5aa44e-d594-11ec-21d6-25d65911031d
using LinearAlgebra,PlutoUI,RDatasets

# ╔═╡ 1f42323f-5074-41f8-bd6a-1806c123e149
md"""
# 5월 17일
"""

# ╔═╡ 32bd60e5-2872-430a-8053-421b6e454f83
# html"""
# <div style="display: flex; justify-content: center;">
# <div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
# <iframe src=
# "
# https://www.youtube.com/embed/playlist?list=PLQqh36zP38-w4vySc_CZBMkKdkXfaxnbV
# "
# width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
# """
	

# ╔═╡ e36c494d-0e97-4010-a57d-fcd60ab93900
md"""
## usings
"""

# ╔═╡ 107c213d-b88e-4138-9a1e-0fba804f6e09
PlutoUI.TableOfContents()

# ╔═╡ ece414e2-e7dd-40f5-b401-5f5039d9b804
md"""
## SVD의 응용: 주성분분석
"""

# ╔═╡ b3ab776c-9dbc-4e17-9f28-0c8f26bba169
md"""
`-` 주성분분석(PCA)은 차원축소에 이용되거나 다중공선성의 해결에 이용되는 기법이다. 
"""

# ╔═╡ 174ef9ba-a3c0-4040-8d31-1cd6962b7b8f
md"""
`-` 우선은 차원축소의 관점에서만 PCA를 서술해보자. 
"""

# ╔═╡ d56fac91-c4c1-4529-9a75-8f4e007cf258
md"""
### PCA의 목적
"""

# ╔═╡ 45ca0b83-4ca4-4137-9cbf-a03ca32d4897
md"""
`-` 데이터: 아래와 같은 $n\times p$ 매트릭스 (혹은 데이터프레임)이 있다고 하자. 

${\bf X}_{n\times p}$

- 통계학과에서 많이 쓰는 그 design matrix라고 생각하면 된다. 
- 여기에서 $n>p$ 를 가정한다. 

"""

# ╔═╡ 762c7cf9-73ab-4fce-9f74-177140560eec
md"""
`-` 소망: (1) ${\bf X}$가 가지고 있는 정보는 거의 유지하면서 (2) 저장비용은 좀 더 저렴한 매트릭스 ${\bf Z}$를 만들고 싶다.
- 정보를 거의 유지한다는 것이 무슨 말? 
- 저장비용이 저렴하다는 것은 무슨 말?

"""

# ╔═╡ 989e6e04-1d01-41be-9d77-973089e0264b
md"""
`-` 소망을 좀 더 구체적으루 수식화하면 아래와 같다. 

(2) ${\bf Z}$의 차원이 $n\times q$ 이며 $q<p$. 

(1) ${\bf Z}$에서 적당한 변환 ${\bf B}$를 하면 언제든지 ${\bf X}$와 비슷한 매트릭스로 복원(reconstruction) 할 수 있다. 즉

${\bf Z}{\bf B} \approx {\bf X}$

을 만족하는 적당한 ${\bf B}$가 존재한다. (이렇게 되면 변환 ${\bf B}$의 치원은 $q \times p$가 되어야 한다.)

"""

# ╔═╡ 9c968667-0a77-42ce-b381-8aa3cf5e31e2
md"""
### SVD를 이용하여 ${\bf Z}$를 찾아보자.
"""

# ╔═╡ 0c36fb86-88e4-4f92-8b34-94944150f273
md"""
`-` 임의의 매트릭스 ${\bf X}_{n\times p}$에 대하여 아래가 성립한다. 단 $n>p$ 를 가정한다. 

${\bf X}_{n\times p}={\bf U}{\bf D}{\bf V}^\top=\sum_{j=1}^{p}U_jd_jV_j^\top$
"""

# ╔═╡ 9c64345f-1dc0-451d-bfcf-13df62a9b525
md"""
`-` 그런데 적당한 $q<p$ 에 대하여 아래가 만족한다고 가정하자. 

${\bf X}_{n \times p}\approx \sum_{j=1}^{q}U_jd_jV_j^\top$
"""

# ╔═╡ 926acde9-9381-4582-85a8-b1d620b3fbf9
md"""
`-` 수식을관찰

$\sum_{j=1}^{q}U_jd_jV_j^\top= 
\begin{bmatrix}
U_1 & U_2 & \dots & U_q 
\end{bmatrix}
\begin{bmatrix}
d_1 & 0 & \dots & 0 \\
0 & d_2 & \dots & 0 \\
0 & 0 & \dots & 0 \\
0 & 0 & \dots & d_q 
\end{bmatrix}
\begin{bmatrix}
V_1^\top \\ V_2^\top \\ \dots \\ V_q^\top
\end{bmatrix}= {\bf \tilde U}{\bf \tilde D}{\bf \tilde V}^\top$
"""

# ╔═╡ b8143b05-9947-4b34-b611-2662ec601ac1
md"""
`-` 따라서 정리하면 아래와 같다. 

${\bf X}_{n \times p}\approx{\bf \tilde U}{\bf \tilde D}{\bf \tilde V}^\top$
"""

# ╔═╡ ffd6d9f2-e4a5-4f06-a75b-242cb4ee71d1
md"""
`-` 여기에서 ${\bf \tilde U}{\bf \tilde D}$의 차원은 $n\times q$가 된다. 이 매트릭스를 ${\bf Z}$로 잡아보자. 즉 

${\bf Z}={\bf \tilde U}{\bf \tilde D}$

이다. 이렇게 되면 ${\bf Z}$의 차원은 $n\times q$ 가 되면서 적당한 변환 ${\bf B}={\bf \tilde V}^\top$를 취하면 ${\bf Z}{\bf B}\approx {\bf X}$가 되므로 PCA의 소망에서 언급한 (1),(2)의 조건을 만족하게 된다. 
"""

# ╔═╡ dd3262c1-d3fb-4feb-ac5b-187833541b19
md"""
### 예제: iris data
"""

# ╔═╡ cee6edd8-a040-4c7c-a417-854904fb413b
md"""
`-` iris 데이터 로드
"""

# ╔═╡ ba4fdd2e-f83b-4b85-9ec2-85119c97ab15
iris = dataset("datasets","iris")

# ╔═╡ 0017fd0e-2ba5-4c30-a2f1-b0f6edf0c3c2
md"""
`-` 각 열에 접근
"""

# ╔═╡ c9d2ece5-6594-454c-9e7a-20a430a08751
iris.SepalLength # 방법1

# ╔═╡ f5acceb2-3523-4f66-ae54-4ea8fb01ee12
iris[:,1] # 방법2

# ╔═╡ 223cf690-d98b-4323-aeb2-43bbf1fc7eb0
md"""
`-` ${\bf X}$ 설정 (이 매트릭스를 차원축소하고 싶다. 거의 손실없이)
"""

# ╔═╡ a94a54f3-48ce-4161-b567-1c46ad0d32dd
X = Array(iris[:,1:4])

# ╔═╡ 238de6fe-735d-4ada-9543-3b78ff3f4b2a
md"""
`-` SVD를 수행
"""

# ╔═╡ 27548de8-5e9a-407f-b3e4-97b596046dfb
U,d,V = svd(X)

# ╔═╡ 28a913da-2deb-4bc1-ba33-4ae91450b98c
md"""
`-` $q=2$로 차원축소 -> $\bf Z$의 dim은 $n\times 2$ 
"""

# ╔═╡ fd4dc6a4-26df-42f0-8fdc-bc9348a0fecf
begin
	Ũ = U[:,1:2]
	D̃ = Diagonal(d[1:2])
	Z = Ũ * D̃
end

# ╔═╡ 4f763d0c-8e57-4bf6-b843-930f23a54565
md"""
`-` reconstruction
"""

# ╔═╡ 0a03f6cd-19aa-4929-b68b-cf17b4fd8744
begin
	B = V[:,1:2]' 
	X̂ = Z * B 
end 

# ╔═╡ a3b010a1-6fa3-426b-8139-47e297d2ea63
md"""
`-` 원본과 비교
"""

# ╔═╡ 48e1efd4-cb2e-4a82-89cb-99b292e9e475
[X X̂]

# ╔═╡ e59a3dca-4055-47bd-9595-5ad28ff37e1c
md"""
### ${\bf Z}$를 계산하는 방법!
"""

# ╔═╡ b7ead0c1-079d-48ee-bf4c-aa699bb50161
md"""
#### 방법 1 (방금 우리가 한 것)
"""

# ╔═╡ 12327076-c21c-4f09-aa3e-28ad52dc8576
md"""
`-` 요약 
"""

# ╔═╡ 23716378-d243-4390-ac47-0d6dbdd09ffd
md"""
(1) ${\bf X}$의 svd를 구한다. 즉 

${\bf X}={\bf U}{\bf D}{\bf V}^\top$ 

를 구한다. 
"""

# ╔═╡ d80fef6b-0cdd-4388-9f43-2e4263746bd4
md"""
(2) ${\bf Z}={\bf \tilde U}{\bf \tilde D}$를 계산한다. 이때 
-  ${\bf \tilde U}=[U_1 \dots U_q]$
-  ${\bf \tilde D}=diag(d_1,\dots,d_q)$
이다 
"""

# ╔═╡ fea795ff-dbcd-4e41-84af-3346167ed6e3
md"""
#### 방법 2
"""

# ╔═╡ f68cb66d-5052-4950-868d-6d56249f0886
md"""
(1) ${\bf X}^\top {\bf X}$를 계산한다.
"""

# ╔═╡ 9f24f9b6-5056-4e0e-843b-726049e9464b
X'X

# ╔═╡ 85bd8e02-96dc-4729-a299-c9ee47210728
md"""
(2) ${\bf X}^\top{\bf X}$의 svd를 구한다. 편의상 아래의 기호로 정의하자. 

${\bf X}^\top{\bf X} ={\bf \Psi} {\bf \Lambda} {\bf \Psi}^\top$
- 이때 왜 ${\bf X}^\top{\bf X} ={\bf \Psi} {\bf \Lambda} {\bf \Psi}^\top$ 와 같이 놓을 수 있지? 
"""

# ╔═╡ 3a016304-9990-4f1a-8974-692f49058376
Ψ,Λ,Ψ = svd(X'X)

# ╔═╡ 32b22a65-445a-4c0a-9f9f-2f916fd6c2ff
md"""
(3) ${\bf Z}={\bf X}{\bf \tilde \Psi}$ 라고 둔다.
- 이때 ${\bf \tilde{\Psi}}=[\Psi_1 ~ \Psi_2]$. 
"""

# ╔═╡ 49825688-e301-4621-a862-5d3e3bc51c96
begin 
	Ψ̃ = Ψ[:,1:2]
	Z′ = X * Ψ̃
end 

# ╔═╡ 68ce9024-463f-4371-a138-962a5251346b
md"""
`-` 확인
"""

# ╔═╡ 6562bfb4-afeb-47f0-928e-d18d4799e924
begin
	B′ = Ψ[:,1:2]' 
	X̂′ = Z′ * B′
	X, X̂, X̂′ 
end 

# ╔═╡ 5f10244a-45a4-451a-8254-b61b77512f96
md"""
## 숙제
"""

# ╔═╡ bfe6ea18-5586-4aa3-966d-28b02d0ca09d
md"""
아이리스 자료의 ${\bf X}$를 (150,3)의 차원을 가지는 ${\bf Z}$로 차원축소하라. 방법1과 방법2를 모두 이용하라. 
""


# ╔═╡ 70f627fa-d3fb-464b-b041-39dab3df6f8f
md"""
## ${\bf X}$의 생성
"""

# ╔═╡ 2d252937-b403-4260-ad4e-b72a7591f201
begin
	μ₁ = [1,2,-2]
	μ₂ = [1,-1,0]
	Σ₁ = [1 0.1 0.5; 0.1 1 0.3 ; 0.5 0.3 1]
	Σ₂ = [1 0.5 0.1; 0.5 1 0.3 ; 0.1 0.3 1]
	dist1= MvNormal(μ₁,Σ₁)
	dist2= MvNormal(μ₂,Σ₂) 
	X = [rand(dist1,10^4) rand(dist2,10^4)]'
end

# ╔═╡ 3eeedf2b-9dae-4294-8c90-08f5d086bdfa
scatter3d(X[:,1],X[:,2],X[:,3],markersize=0.2)

# ╔═╡ 957dbbb5-fee6-4346-980f-70aa55f44e90
md"""
## ${\bf Z}$의 계산
"""

# ╔═╡ c3263c6d-0e1c-4bbf-85b2-b1a3e155496c
md"""
`-` 목표: ${\bf X}_{n\times 3}$를 차원축소하여 ${\bf Z}_{n\times 2}$로 만들자. (각 observation을 3차원에서 2차원으로 변환)
"""

# ╔═╡ 1c47fdea-a511-43e6-9a2e-bdd044560682
md"""
### 방법1: 직접계산 (${\bf Z}={\bf \tilde U}{\bf \tilde D}$)
"""

# ╔═╡ 83e16950-3d5f-4eeb-a15a-46855a16ddc6
let
	U,d,V = svd(X)
	Ũ = U[:,1:2] # p=3, q=2 설정 
	D̃ = Diagonal(d[1:2])
	Z=Ũ*D̃ 
end

# ╔═╡ 4057a796-d001-4a91-b877-adfd7240b903
md"""
### 방법2: ${\bf Z}={\bf X}{\bf \tilde V}$ 이용
"""

# ╔═╡ 6ed3b3e5-6add-40ca-b58c-dc7a53389c6b
md"""
`-` 관찰: ${\bf Z}={\bf \tilde U}{\bf \tilde D}={\bf X}{\bf \tilde V}$
- 증명: [보충자료(pdf)](https://github.com/guebin/SC2022/blob/main/SC2022_0519_supp.pdf)
"""

# ╔═╡ e4251a72-bab1-43ac-8cc3-08d7e24405b9
let 
	U,d,V = svd(X)
	Ṽ = V[:,1:2] 
	X*Ṽ
end 

# ╔═╡ ac09e9fe-9e6e-4d29-b3a6-70296f2ec825
md"""
`-` 결론: ${\bf Z}={\bf \tilde U}{\bf \tilde D}={\bf X}{\bf \tilde V}$ 이므로, ${\bf \tilde U}{\bf \tilde D}$를 계산해서 ${\bf Z}$를 얻을 수도 있고 ${\bf X}{\bf \tilde V}$를 계산하여 ${\bf Z}$를 얻을 수도 있다.
"""

# ╔═╡ 137f5a6b-2c47-4359-9a30-565a1575e9bd
md"""
## ${\bf X}={\bf U}{\bf D}{\bf V}^\top$에서 ${\bf V}$의 계산방법
"""

# ╔═╡ 84ae5e8f-39e5-42fd-b5f3-7ca2b3009702
md"""
### 방법1: 직접계산
"""

# ╔═╡ 6901a2a9-163e-4602-a7a1-58fb6d2c1e44
let
	U,d,V = svd(X)
	V
end

# ╔═╡ 321d4559-022f-4ec4-8ef9-288af6973803
md"""
### 방법2: ${\bf X}'{\bf X}$ 고유값분해 이용
"""

# ╔═╡ 60a9d85f-6117-4540-b1d0-359b8cf0ccff
let
	λ,Ψ = eigen(X'X)
end

# ╔═╡ 0f200543-8a25-4b7d-92f4-c5972eccaa96
md"""
- 결과를 보니까 ${\bf V}$와 ${\bf \Psi}$는 column을 나열한 순서가 조금 다르고 (하나는 대각행렬에 해당하는 원소를 큰순서대로 나열하고 하나는 작은 순서대로 나열해서 그렇다) column의 부호도 조금 다르게 나오긴 한다. 그래서 ${\bf V}$와 ${\bf \Psi}$가 완전히 같은것은 아니다. 
- 하지만 이러한 사소한 차이는 별로 상관이 없다. 순서가 다른것은 맞춰주면 되고 부호는 달라져도 상관없다. 
"""

# ╔═╡ ec1262be-0742-4b09-9543-5888f5f695c3
md"""
### 방법3: ${\bf X}'{\bf X}$의 SVD를 이용 
"""

# ╔═╡ 3e559a77-5869-49bf-b190-98e3557a9a79
let 
	UU,dd,VV = svd(X'X)
	VV
end 

# ╔═╡ 3e572717-4415-4926-831f-88fb645136d6
let 
	UU,dd,VV = svd(X'X)
	UU
end 

# ╔═╡ 22c14226-96bc-4610-ba22-cf6e9bc9e684
md"""
## 정리: ${\bf X} \to {\bf Z}  \to {\bf \hat X}$ 의 다양한 구현방법
"""

# ╔═╡ 70c85338-d303-4b1c-8822-703730639e1d
md"""
### 코드1: ${\bf Z}={\bf \tilde U}{\bf \tilde D} \to {\bf \hat X}={\bf Z}{\bf \tilde V}^\top$ (X의 SVD를 이용)
"""

# ╔═╡ 27b6f9c6-47b8-412c-8e16-0565d9c2bfb2
let
	U, d, V = svd(X)
	Z = U[:,1:2] * Diagonal(d[1:2])
	X̂ = Z*V[:,1:2]'
	[X X̂]
end

# ╔═╡ 56fdfb52-42a3-4d60-8e4f-e799b9f250db
md"""
### 코드2: ${\bf Z}={\bf X}{\bf \tilde V} \to {\bf \hat X}={\bf Z}{\bf \tilde V}^\top$ (X의 SVD를 이용)
"""

# ╔═╡ 8199694b-099c-48a2-bfa5-7b2d7850bf93
let
	U, d, V = svd(X)
	Z = X * V[:,1:2]
	X̂ = Z*V[:,1:2]'
	[X X̂]
end

# ╔═╡ 71bcc0e1-258f-4c68-8c96-a8a9e891a5a8
md"""
### 코드3: ${\bf Z}={\bf X}{\bf \tilde V} \to {\bf \hat X}={\bf Z}{\bf \tilde V}^\top$ (X'X의 고유값분해를 이용)
"""

# ╔═╡ 94ff4be3-45f0-4c8d-beb1-9272c89732c2
let
	λ, Ψ = eigen(X'X)
	V = [Ψ[:,3] Ψ[:,2] Ψ[:,1]]
	Z = X * V[:,1:2]
	X̂ = Z * V[:,1:2]'
	[X X̂]
end

# ╔═╡ 4b6032ac-02b2-459f-a253-eee7e0627a29
md"""
### 코드4: ${\bf Z}={\bf X}{\bf \tilde V} \to {\bf \hat X}={\bf Z}{\bf \tilde V}^\top$ (X'X의 SVD를 이용)
"""

# ╔═╡ 08fb4c56-787f-4f63-8301-60946682617a
let
	UU,dd,VV = svd(X'X)
	Z = X * VV[:,1:2]
	X̂ = Z * VV[:,1:2]'
	[X X̂]
end

# ╔═╡ d7dbbc18-ca88-4688-8ebf-645df5227e94
let
	UU,dd,VV = svd(X'X)
	Z = X * UU[:,1:2]
	X̂ = Z * UU[:,1:2]'
	[X X̂]
end

# ╔═╡ 46c23937-edbb-4a41-a5f5-158496384524
md"""
## ${\bf Z}$를 잘만들었다는 의미?
"""

# ╔═╡ a0a2c132-db69-48f2-98a3-b5c21f5a6cbe
scatter3d(X[:,1],X[:,2],X[:,3],markersize=0.2)

# ╔═╡ fcf3fcfb-d40e-4871-ad22-07d81889dd38
md"""
-  ${\bf X}$에서 데이터를 관찰하면 점들이 두개의 군집으로 이루어졌다는 "정보"를 확인할 수 있음. 
- 이 정보는 ${\bf Z}$에서 관찰해도 사라지지 않아야 함. 
"""

# ╔═╡ 9f71e05e-aed0-4819-a70a-575b655896b1
md"""
-  ${\bf Z}$에서도 두개의 군집에 대한 정보는 잘 유지되고 있다.
"""

# ╔═╡ 292cc2d2-d98f-4ed3-a2f4-6ff180fc5388
md"""
## V의 column에서 부호를 마음대로 변경
"""

# ╔═╡ 4268bb5f-848e-469a-bdc3-74485b1b7507
md"""
`-` V1의 column의 부호를 임의로 바꿔서 넣으면 당연히 ${\bf Z}$의 값 자체는 달라진다. 
"""

# ╔═╡ 674582c7-cf22-4de8-b28f-37f0efd2c8d1
begin 
	UU,dd,VV = svd(X'X) 
	_VV = [VV[:,1] -VV[:,2] VV[:,3]] # 두번쨰 컬럼의 부호를 임의로 바꿔버렸음 
	Z = X * VV[:,1:2] # 부호를 안바꾼 V1으로 Z를 구함 
	_Z = X * _VV[:,1:2] # 부호를 바꾼 V2으로 Z를 구함 
	[Z _Z]	
end 

# ╔═╡ 2669bec0-b269-4132-a9f9-19017f820b35
scatter(Z[:,1],Z[:,2],markersize=0.2)

# ╔═╡ 14f043b3-9237-45b4-949a-d9f13d6b382a
begin
	X̂ = Z * VV[:,1:2]'
	_X̂ = _Z * _VV[:,1:2]' 
	[X̂ _X̂] # 그런데 reconstruction은 같다. (왜?) <-- 다음에 풀게요 
end 

# ╔═╡ 87675737-6a0b-4358-bd33-0e3e329da94f
md"""
`-` 결론: 이 예제에서 Z와 Z′는 그 값이 똑같진 않지만 
-  ${\bf X}$가 가지고 있는 정보는 그대로 유지하면서 
- 저장비용은 좀더 저렴한 매트릭스를 만들고 싶다. 
라는 관점에서는 같은 매트릭스로 봐도 된다. 
"""

# ╔═╡ 9568ac40-783a-4b8b-ad19-c1bf0952d4ad
md"""
## 숙제
"""

# ╔═╡ 33ed6d3e-26da-4b19-bcc2-5f22135eb1d4
md"""
`-` `scatter(_Z[:,1],_Z[:,2],markersize=0.2)`를 이용하여 스캐터플랏을 그려라. 두 개의 그룹이 여전히 잘 구분되는지 확인하라.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
RDatasets = "ce6b1742-4840-55fa-b093-852dadbb1d8b"

[compat]
PlutoUI = "~0.7.58"
RDatasets = "~0.7.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "a19025e7c9ef3323cc1d57e271b64b306fa491f9"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "297b6b41b66ac7cbbebb4a740844310db9fd7b8c"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "6c834533dc1fabd820c1db03c839bf97e45a3fab"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.14"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "1568b28f91293458345dabba6a5ea3f183250a61"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.8"

    [deps.CategoricalArrays.extensions]
    CategoricalArraysJSONExt = "JSON"
    CategoricalArraysRecipesBaseExt = "RecipesBase"
    CategoricalArraysSentinelArraysExt = "SentinelArrays"
    CategoricalArraysStructTypesExt = "StructTypes"

    [deps.CategoricalArrays.weakdeps]
    JSON = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
    SentinelArrays = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
    StructTypes = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "97d79461925cdb635ee32116978fc735b9463a39"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.19"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExprTools]]
git-tree-sha1 = "27415f162e6028e81c72b82ef756bf321213b6ec"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.10"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "82d8afa92ecf4b52d78d869f038ebfb881267322"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "9f00e42f8d99fdde64d40c8ea5d14269a2e2c1aa"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.21"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

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

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

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
version = "2.28.2+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "4cc0c5a83933648b615c36c2b956d94fda70641e"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.7"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

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

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "88b895d13d53b5577fd53379d913b9ab9ac82660"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.RData]]
deps = ["CategoricalArrays", "CodecZlib", "DataFrames", "Dates", "FileIO", "Requires", "TimeZones", "Unicode"]
git-tree-sha1 = "19e47a495dfb7240eb44dc6971d660f7e4244a72"
uuid = "df47a6cb-8c03-5eed-afd8-b6050d6c41da"
version = "0.8.3"

[[deps.RDatasets]]
deps = ["CSV", "CodecZlib", "DataFrames", "FileIO", "Printf", "RData", "Reexport"]
git-tree-sha1 = "2720e6f6afb3e562ccb70a6b62f8f308ff810333"
uuid = "ce6b1742-4840-55fa-b093-852dadbb1d8b"
version = "0.7.7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "0e7508ff27ba32f26cd459474ca2ede1bc10991f"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TZJData]]
deps = ["Artifacts"]
git-tree-sha1 = "b69f8338df046774bd975b13be9d297eca5efacb"
uuid = "dc5dba14-91b3-4cab-a142-028a31da12f7"
version = "1.1.0+2023d"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TimeZones]]
deps = ["Artifacts", "Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "Scratch", "TZJData", "Unicode", "p7zip_jll"]
git-tree-sha1 = "cc54d5c9803309474014a8955a96e4adcd11bcf4"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.14.0"

    [deps.TimeZones.extensions]
    TimeZonesRecipesBaseExt = "RecipesBase"

    [deps.TimeZones.weakdeps]
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.TranscodingStreams]]
git-tree-sha1 = "71509f04d045ec714c4748c785a59045c3736349"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.7"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

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

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─1f42323f-5074-41f8-bd6a-1806c123e149
# ╠═32bd60e5-2872-430a-8053-421b6e454f83
# ╟─e36c494d-0e97-4010-a57d-fcd60ab93900
# ╠═3c5aa44e-d594-11ec-21d6-25d65911031d
# ╠═107c213d-b88e-4138-9a1e-0fba804f6e09
# ╟─ece414e2-e7dd-40f5-b401-5f5039d9b804
# ╟─b3ab776c-9dbc-4e17-9f28-0c8f26bba169
# ╟─174ef9ba-a3c0-4040-8d31-1cd6962b7b8f
# ╟─d56fac91-c4c1-4529-9a75-8f4e007cf258
# ╟─45ca0b83-4ca4-4137-9cbf-a03ca32d4897
# ╟─762c7cf9-73ab-4fce-9f74-177140560eec
# ╟─989e6e04-1d01-41be-9d77-973089e0264b
# ╟─9c968667-0a77-42ce-b381-8aa3cf5e31e2
# ╟─0c36fb86-88e4-4f92-8b34-94944150f273
# ╟─9c64345f-1dc0-451d-bfcf-13df62a9b525
# ╟─926acde9-9381-4582-85a8-b1d620b3fbf9
# ╟─b8143b05-9947-4b34-b611-2662ec601ac1
# ╟─ffd6d9f2-e4a5-4f06-a75b-242cb4ee71d1
# ╟─dd3262c1-d3fb-4feb-ac5b-187833541b19
# ╟─cee6edd8-a040-4c7c-a417-854904fb413b
# ╠═ba4fdd2e-f83b-4b85-9ec2-85119c97ab15
# ╟─0017fd0e-2ba5-4c30-a2f1-b0f6edf0c3c2
# ╠═c9d2ece5-6594-454c-9e7a-20a430a08751
# ╠═f5acceb2-3523-4f66-ae54-4ea8fb01ee12
# ╟─223cf690-d98b-4323-aeb2-43bbf1fc7eb0
# ╠═a94a54f3-48ce-4161-b567-1c46ad0d32dd
# ╟─238de6fe-735d-4ada-9543-3b78ff3f4b2a
# ╠═27548de8-5e9a-407f-b3e4-97b596046dfb
# ╟─28a913da-2deb-4bc1-ba33-4ae91450b98c
# ╠═fd4dc6a4-26df-42f0-8fdc-bc9348a0fecf
# ╟─4f763d0c-8e57-4bf6-b843-930f23a54565
# ╠═0a03f6cd-19aa-4929-b68b-cf17b4fd8744
# ╟─a3b010a1-6fa3-426b-8139-47e297d2ea63
# ╠═48e1efd4-cb2e-4a82-89cb-99b292e9e475
# ╟─e59a3dca-4055-47bd-9595-5ad28ff37e1c
# ╟─b7ead0c1-079d-48ee-bf4c-aa699bb50161
# ╟─12327076-c21c-4f09-aa3e-28ad52dc8576
# ╟─23716378-d243-4390-ac47-0d6dbdd09ffd
# ╟─d80fef6b-0cdd-4388-9f43-2e4263746bd4
# ╟─fea795ff-dbcd-4e41-84af-3346167ed6e3
# ╟─f68cb66d-5052-4950-868d-6d56249f0886
# ╠═9f24f9b6-5056-4e0e-843b-726049e9464b
# ╟─85bd8e02-96dc-4729-a299-c9ee47210728
# ╠═3a016304-9990-4f1a-8974-692f49058376
# ╟─32b22a65-445a-4c0a-9f9f-2f916fd6c2ff
# ╠═49825688-e301-4621-a862-5d3e3bc51c96
# ╟─68ce9024-463f-4371-a138-962a5251346b
# ╠═6562bfb4-afeb-47f0-928e-d18d4799e924
# ╟─5f10244a-45a4-451a-8254-b61b77512f96
# ╟─bfe6ea18-5586-4aa3-966d-28b02d0ca09d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
