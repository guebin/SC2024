### A Pluto.jl notebook ###
# v0.19.42

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

# ╔═╡ 128613c3-baba-426a-adbf-8abed179eb49
using PlutoUI,Plots,HTTP,CSV,DataFrames,LinearAlgebra,Statistics,Random

# ╔═╡ 7ec1510e-130c-11ef-1e0e-518d9cc440ea
md"""
# 12wk-1: 능형회귀 (1)
"""

# ╔═╡ a366d78e-5826-4a65-a3fe-3469fee80f9a
md"""
## 1. 강의영상
"""

# ╔═╡ 4c8cb126-08bd-478e-8472-fd1547e0d128
html"""
<div style="display: flex; justify-content: center;">
<div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
<iframe src=
"
https://youtube.com/embed/playlist?list=PLQqh36zP38-y5K74hZNfuDUv2WsUHehw_&si=NO3AvV5GHNHRh65c
"
width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ 2d2f8c88-07f9-4e8b-86ae-fc9828aba187
md"""
## 2. Imports
"""

# ╔═╡ 3317fffc-955d-4073-9123-17688a4c6f61
PlutoUI.TableOfContents()

# ╔═╡ 18996e8e-c802-4e03-a3c6-c4d83aa5afdc
Plots.plotly()

# ╔═╡ 31f2972f-c83f-42ce-999c-842c157899d3
md"""
## 3. 다중공선성
"""

# ╔═╡ b5277d9e-216f-444f-bcae-e7ddf5d2e9b5
df = DataFrame(CSV.File(HTTP.get("https://raw.githubusercontent.com/guebin/SC2024/main/toeic.csv").body))

# ╔═╡ df19fb9c-8f86-4d55-9899-c868c3512411
n = 5000

# ╔═╡ f9f88913-a23a-4f02-a1cd-f255feb8c834
begin 
	X1,X2,X3 = eachcol(df)
	X = [X1 X2 X3] 
end

# ╔═╡ b45868ea-681e-41e3-b9c0-2d0c918e574d
begin 
	Random.seed!(43052) 
	β1 = 600 
	β2 = 5 
	β3 = 0
	β = [β1, β2, β3]
	σ = 300
	ϵ = σ*randn(n)
	#y = β1*X1 + β2*X2 + β3*X3 + ϵ
	y = X*β + ϵ
end

# ╔═╡ d055be6a-aef8-47d0-b0a2-a56b85089bfc
md"""
## 4. 능형회귀
"""

# ╔═╡ 8b6bafd7-7dd9-4c2b-a059-d18fc2ad22f8
md"""
### A. 문제정리
"""

# ╔═╡ 57b42243-1fa1-4a69-b116-044b30d7819b
md"""
-- 일반적인 회귀문제는 아래를 푸는 것이다. 
"""

# ╔═╡ b70e73b7-c2ae-4c94-9f85-ad0195b5858a
md"""
$$\underset{{\boldsymbol \beta} \in \mathbb{R}^p}{\operatorname{argmin}}\bigg\{ \big({\bf y}-{\bf X}{\boldsymbol \beta} \big)^\top \big({\bf y}-{\bf X}{\boldsymbol \beta}  \big)\bigg\}$$
"""

# ╔═╡ 6601a37b-cb7c-417d-a438-c5c18c19e0cf
md"""
-- 이에 대한 수학적인 해는 $\hat{\boldsymbol \beta}=\big({\bf X}^\top {\bf X}\big)^{-1}{\bf X}^\top{\bf y}$ 임을 너무나도 잘 알고 있지만, 우리의 예제에서는 이 수학적인 해가 별로 쓸모가 없다는 사실을 확인했다. 왜 이런일이 생길까?
"""

# ╔═╡ 23d619d0-e790-4ddb-97ec-735e42fa85bc
md"""
-- 편의상 GPA에 대한 추정값 600은 정확하게 추정했다고 가정하자. 즉 이제 아래의 모형을 가정한다.

$\tilde{\bf y}= {\bf y}-\beta_1{\boldsymbol X}_1 =\beta_2{\boldsymbol X}_2 + \beta_3{\boldsymbol X}_3 +{\boldsymbol \epsilon}$
"""

# ╔═╡ 68d0574b-dbd6-4b68-920e-e672fe876ec6
ỹ = y - β1*X1

# ╔═╡ ac6df6a5-6742-446a-b119-a0a685de2f3f
md"""
-- 손실함수를 정의하자.
"""

# ╔═╡ cc6f7ff7-e7ac-409e-9f63-145c58e4ade5
loss(β2,β3) = (ỹ-β2*X2-β3*X3)'*(ỹ-β2*X2-β3*X3)/n

# ╔═╡ f05be255-fae4-4421-a920-c863d0c389a6
md"""
-- 손실함수를 그려보자.
"""

# ╔═╡ cbb62ff1-fe39-468e-91a6-0431192fb25d
begin
	β̂2s = -10:0.5:15
	β̂3s = -10:0.5:15
	p1 = plot(β̂2s,β̂3s,loss,st=:surface,colorbar=false,alpha=0.9)
	p2 = plot(β̂2s,β̂3s,loss,st=:contour,colorbar=false,levels=100)
	plot(p1,p2)
end

# ╔═╡ 3bad3692-7f76-494c-93d5-f3951331c74e
md"""
- ``\hat{\beta}_2 + \hat{\beta}_3 \approx 5`` 이라면 loss값을 거의 최소값과 비슷하게 만든다. 
"""

# ╔═╡ 758668d4-7d98-46db-8f8a-9a81abe6ad25
md"""
-- 수식적으로 그럴듯해 보여도, 모두 바람직한 모형은 아니다. 아래는 모두 참모형이라 생각되어지는 상황이다. 

1. ``\hat{\beta}_2=2.5,\quad \hat{\beta}_3=2.5.``
2. ``\hat{\beta}_2=5,\quad \hat{\beta}_3=0.``
3. ``\hat{\beta}_2=100,\quad \hat{\beta}_3=-95.``

그렇지만 상식적으로 3은 용납할 수 없다. 계수값은 최소한 0~5 사이의 값이었으면 좋겠다. (-0.1정도는 괜찮을 듯)
"""

# ╔═╡ a82b89e1-ac4b-45e0-9e46-db58bc303812
md"""
-- 각 경우에 대하여 loss를 계산해보자.
"""

# ╔═╡ a8e9091b-aa4b-49c8-804f-b0db309329a8
let
	@show loss(2.5,2.5)
	@show loss(5,0)
	@show loss(100,-95)
end 

# ╔═╡ cbe878f4-07b2-4f21-8a63-7468b2633b6b
md"""
- 이대로라면 컴퓨터입장에서는 경우3이 가장 좋다고 생각하겠는걸?
"""

# ╔═╡ 90ff8925-e0d6-4982-b19d-0f78b5016545
md"""
### B. 해결책
"""

# ╔═╡ 85ad9c84-87c8-49d9-9139-63605442aa2d
md"""
!!! warning "아이디어: 추정된 계수의 절대값이 크면 벌점을 주자"
	아래와 같은 계수추정값들이 있다고 하자.

	1. ``(\hat{\beta}_2,\hat{\beta}_3) = (0,5)``
	2. ``(\hat{\beta}_2,\hat{\beta}_3) = (1,4)``
	3. ``(\hat{\beta}_2,\hat{\beta}_3) = (2,3)``
	4. ``(\hat{\beta}_2,\hat{\beta}_3) = (2.5,2.5)``
	5. ``(\hat{\beta}_2,\hat{\beta}_3) = (3,2)``
	6. ``(\hat{\beta}_2,\hat{\beta}_3) = (4,1)``
	7. ``(\hat{\beta}_2,\hat{\beta}_3) = (5,0)`` -- 찐 참모형

	이러한 추정값들은 대충 합리적이라 여겨지며 우리는 이러한 값을 사실 추정하고 싶다. (7이 진짜 참모형이지만 1-6도 나쁘지 않은 선택임) 그럼 이제 아래의 값들을 살펴보자. 

	8. ``(\hat{\beta}_2,\hat{\beta}_3) = (-5,10)``
	9. ``(\hat{\beta}_2,\hat{\beta}_3) = (-95,100)``
	10. ``(\hat{\beta}_2,\hat{\beta}_3) = (-995,1000)``
	11. ``(\hat{\beta}_2,\hat{\beta}_3) = (-9995,10000)``

	이 값들은 딱 봐도 짜증나는 상황(그렇지만 수학적으로 가능할 것 같은 상황)이며, 우리는 이러한 값들 추정하고 싶지 않다. 우리가 원하는 추정값과 원하지 않는 추정값을 관찰하면서 얻은 직관은 아래와 같다. 

	> ``\hat{\beta}_2`` 와 ``\hat{\beta}_3`` 의 절대값이 커질수록 우리가 원하는 추정량은 아닌것 같다.

	**그렇다면 ``\hat{\beta}_2`` 와 ``\hat{\beta}_3`` 의 절대값이 클수록 손실함수를 더 크게 만들어 버리면 어떨까??!**

"""

# ╔═╡ 6216c9c1-3d3a-4e2c-9fd7-d18cba3b6147
md"""
-- ``\lambda \geq 0`` 에 대하여 아래와 같은 손실함수를 고려해보자.

$loss_{\text{L}^2} := \big({\bf y}-{\bf X}{\boldsymbol \beta} \big)^\top \big({\bf y}-{\bf X}{\boldsymbol \beta}  \big) + \lambda {\boldsymbol \beta}^\top{\boldsymbol \beta}$
"""

# ╔═╡ 1d20583d-bddb-4c38-9da5-03d9e4c53f03
md"""
-- $L_2$ 벌점을 추가할 경우 아래를 계산해보자.

1. ``\hat{\beta}_2=2.5,\quad \hat{\beta}_3=2.5.``
2. ``\hat{\beta}_2=5,\quad \hat{\beta}_3=0.``
3. ``\hat{\beta}_2=100,\quad \hat{\beta}_3=-95.``

"""

# ╔═╡ d3416f40-9387-4a44-986f-8657c8ea1c00
let
	@show λ
	println("---")
	@show loss(2.5,2.5)
	@show l2(2.5,2.5)
	@show loss_l2(2.5,2.5) # 이제는 이게 제일 작음
	println("---")
	@show loss(0,5)
	@show l2(0,5)
	@show loss_l2(0,5)
	println("---")
	@show loss(100,-95)
	@show l2(100,-95)
	@show loss_l2(100,-95)
end 

# ╔═╡ 324c7de6-c404-470a-afe3-cc756cf9fa94
md"""
### C. 변형된 손실함수 시각화
"""

# ╔═╡ 1109294c-c93c-4156-a301-e76260f5b1a4
md"""
λ = $(@bind λ Slider([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7],show_value=true, default=1e5))
"""

# ╔═╡ 5ba39dcd-1f21-4c83-a420-59de892a8121
begin 
	l2(β̂2,β̂3) = λ*(β̂2^2 + β̂3^2)
	loss_l2(β̂2,β̂3) = loss(β̂2,β̂3) + l2(β̂2,β̂3)
end 

# ╔═╡ 05665ffd-fba5-49c4-bae7-65c6b02ebc64
md"""
*Figure 1*:  $(x,y,z) = \big(\hat{\beta}_2,~ \hat{\beta}_3,~ loss(\hat{\beta}_2,\hat{\beta}_3)\big)$
"""

# ╔═╡ 278ef454-7fcc-4155-89f4-6fd93d168feb
let
	p1 = plot(β̂2s,β̂3s,loss,st=:surface,colorbar=false,alpha=0.9)
	p2 = plot(β̂2s,β̂3s,loss,st=:contour,colorbar=false,levels=100)
	plot(p1,p2)
end 

# ╔═╡ 14f7efbc-448d-44e1-b107-d52f05bc8520
md"""
*Figure 2*:  $(x,y,z) = \Big(\hat{\beta}_2, ~\hat{\beta}_3,~ \lambda(\hat{\beta}_2^2 + \hat{\beta}_3^2)\Big)$
"""

# ╔═╡ c5bdf20e-e89f-46e6-9a9a-2c16a19a54c4
let
	p3 = plot(β̂2s,β̂3s,l2,st=:surface,colorbar=false,alpha=0.9)
	p4 = plot(β̂2s,β̂3s,l2,st=:contour,colorbar=false,levels=100)
	plot(p3,p4)
end

# ╔═╡ 21a75194-76de-43bd-abbe-89d3e28e6d8e
md"""
*Figure 3*:  $(x,y,z) = \Big(\hat{\beta}_2,~ \hat{\beta}_3,~ loss(\hat{\beta}_2,\hat{\beta}_3) +\lambda(\hat{\beta}_2^2 + \hat{\beta}_3^2)\Big)$
"""

# ╔═╡ 14e02673-0863-4aa3-9ba6-af0bc909d0cc
let
	p5 = plot(β̂2s,β̂3s,loss_l2,st=:surface,colorbar=false,alpha=0.9)
	p6 = plot(β̂2s,β̂3s,loss_l2,st=:contour,colorbar=false,levels=100)
	plot(p5,p6)
end 

# ╔═╡ 6d545544-3516-4560-8887-e4413f3c0fc0
md"""
### D. 해를 구하는 방법
"""

# ╔═╡ f5b49064-9f8b-4238-97a4-780b4023a641
md"""
-- 다중공선성이 의심되는 경우에는 단순하게 아래를 최소화 하는 것 보다

$loss := \big({\bf y}-{\bf X}{\boldsymbol \beta} \big)^\top \big({\bf y}-{\bf X}{\boldsymbol \beta}  \big)$

``\lambda \geq 0`` 에 대하여 아래와 같은 손실함수를 최소화하는 $\hat{\boldsymbol \beta}$ 을 구하는 것이 더 이득인 것 같다. 

$loss_{\text{L}^2} := \big({\bf y}-{\bf X}{\boldsymbol \beta} \big)^\top \big({\bf y}-{\bf X}{\boldsymbol \beta}  \big) + \lambda {\boldsymbol \beta}^\top{\boldsymbol \beta}$
"""

# ╔═╡ 8a47abe4-22de-48fc-b98c-c4cce53502bb
md"""
-- ``loss_{\text{L}^2}``를 최소화하는 수학적인 해는 아래와 같다.
"""

# ╔═╡ b3d527a8-3151-4f30-8ccd-b53b298855de
md"""
$\hat{\boldsymbol \beta}=({\bf X}^\top {\bf X}+\lambda {\bf I})^{-1}{\bf X}^\top {\bf y}$
"""

# ╔═╡ 2fed1588-f9cc-4c6d-a8ff-751352a1fcd0
let
	λ = 50
	β̂ = inv(X'X + λ*I)X'y
end 

# ╔═╡ 61b295ca-be2d-4326-bff5-0fd860d31919
md"""
-- 결과는 그럭저럭 괜찮음. 
- ``\hat{\beta}_3``, ``\hat{\beta}_3`` 에 대한 추정값이 괜찮게 나온것은 긍정적임. 
- 그런데 ``\hat{\beta}_1``의 추정값은 $\lambda$ 값을 키울수록 600보다 작게 추정된다. 
"""

# ╔═╡ 057f12ba-1e54-481a-b9f4-a6e2bcad6805
md"""
-- 다른 평행세계에 대하여서도 ``\hat{\beta}_1``,``\hat{\beta}_2``, ``\hat{\beta}_3`` 를 각각 구해보고 그 분포를 살펴보자. 
"""

# ╔═╡ 392d159d-c7a4-47f8-bdfd-c28fdfc1e6a6
let 
	N = 10000
	E = 300*randn(n,N)
	Y = (600*X1 + 5*X2) .+ E
	λ = 50
	B̂ = inv(X'X + λ*I)X'Y
	β̂1s,β̂2s,β̂3s = eachrow(B̂)
	p1 = histogram(β̂1s,alpha=0.5,label="β̂₁")
	p2 = histogram(β̂2s,alpha=0.5,label="β̂₂")
	p3 = histogram(β̂3s,alpha=0.5,label="β̂₃")
	plot(p1,p2,p3)
end

# ╔═╡ e1fea1fb-7a22-41de-b84a-b02c77c728c9
md"""
- 직관1: $\lambda$를 키울수록 추정량의 분산은 줄어든다. (이건 좋은거)
- 직관2: $\lambda$를 키울수록 추정량이 작은값을 가진다. (이건 나쁜거)
"""

# ╔═╡ 1c21d098-4764-469a-9a93-9c091304b7cd
md"""
## 5. 능형회귀의 특징
"""

# ╔═╡ 42cd97a0-c604-4656-9d4c-fc4a8f268580
md"""
### A. 추정량의 평균
"""

# ╔═╡ cc218bfe-8862-4e48-8302-88b6a17da2a2
md"""
능형회귀로 얻은 추정량 $\hat{\boldsymbol \beta}=({\bf X}^\top{\bf X}+\lambda {\bf I})^{-1}{\bf X}^\top {\bf y}$ 에 대하여 $\mathbb{E}(\hat{\boldsymbol \beta})$의 값을 조사해보자. 
"""

# ╔═╡ fbeb2602-7d61-4a6e-9ca6-bf8618721d1d
md"""
-- 방법1: 이론적으로 조사하자.
"""

# ╔═╡ d7206a25-b7f6-4c59-8efc-268f3485fe2d
md"""
$\mathbb{E}(\hat{\boldsymbol\beta}) =({\bf X}^\top{\bf X}+\lambda {\bf I})^{-1}{\bf X}^\top{\bf X}{\boldsymbol \beta}$
"""

# ╔═╡ ebb19d9a-0ddc-40a5-9533-79639bba2e82
md"""
이 값은 우리가 기대하는 좋은 추정량의 성질인 $\mathbb{E}(\hat{\boldsymbol \beta}) = {\boldsymbol \beta}$ 와 다르다. 
"""

# ╔═╡ b71607cd-dfb5-4138-9042-fdd1b7b93bba
md"""
!!! warning ""
	아쉬움.. ``\lambda =0`` 이었다면 $\mathbb{E}(\hat{\boldsymbol \beta})={\boldsymbol \beta}$ 이었을텐데.. 
"""

# ╔═╡ 096103d0-b7e3-4f60-b7d9-a69a58b4d14e
inv(X'X+0.000000000001*I)X'X*[600,5,0]

# ╔═╡ 0401b893-61ed-426a-8f0c-cfa365c56022
md"""
*Figure: $\lambda$에 따른 추정량의 평균변화 (이론)*
"""

# ╔═╡ 80b71ff3-3b35-4a6d-bacb-5a69d735dabc
let 
	Eβ̂(λ) = inv(X'X+λ*I)X'X*β
	λs = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11]
	plot(λ -> Eβ̂(λ)[1], λs, xscale= :log10, yscale= :log10, label="E(β̂₁)")
	plot!(λ -> Eβ̂(λ)[2], label="E(β̂₂)")
	plot!(λ -> Eβ̂(λ)[3], label="E(β̂₃)")
end

# ╔═╡ 86ee8e65-1220-4d8f-8801-a51e1924cd85
md"""
-- 방법2: 시뮬레이션으로 조사하자.
"""

# ╔═╡ 828c4d10-263b-46d2-b36f-28f6a0f3c612
let 
	λ = 10
	N = 1000
	E = σ*randn(n,N)
	Y = X*β .+ E
	B̂ = inv(X'X+λ*I)X'Y
	j = ones(N)
	B̂ * j/N
end

# ╔═╡ b1266c8b-2be6-4fb0-a791-38d9f2844d0b
md"""
- ``\lambda`` 키울수록 true 값과 멀어진다. 
- ``\lambda`` 키울수록 $\mathbb{E}(\hat{\beta})$의 관점에서는 손해!
"""

# ╔═╡ 101cd18a-0b4b-4e1e-8dee-6d687ebedc30
md"""
*Figure: $\lambda$에 따른 추정량의 평균변화 (시뮬)*
"""

# ╔═╡ c029a13f-ab20-47f0-92e0-edd767302f85
let 
	function Êβ̂(λ)
		N = 1000
		E = randn(n,N)*σ
		Y = X*β .+ E
		B̂ = inv(X'X+λ*I)X'Y
		j = ones(N)
		return B̂ * j/N
	end
	λs = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11]
	plot(λ -> Êβ̂(λ)[1],λs,xscale= :log10,yscale= :log10, label="Ê(β̂₁)")
	plot!(λ -> Êβ̂(λ)[2], label="Ê(β̂₂)")
	plot!(λ -> Êβ̂(λ)[3], label="Ê(β̂₃)")
end

# ╔═╡ 0b7a7efb-e3e2-42f5-840c-e8245660e840
md"""
!!! warning "능형회귀 추정량에 대한 수식어"
	능형회귀로 얻은 추정량 $\hat{\boldsymbol \beta}=({\bf X}^\top{\bf X}+\lambda {\bf I})^{-1}{\bf X}^\top {\bf y}$ 는 *biased estimator* / *shrinkage estimator* 라고 불린다. 이는 모두 $\mathbb{E}(\hat{\boldsymbol \beta})$ 에 대한 성질때문에 생긴 수식어이다. 여기에서 *biased* 라는 의미는 $\mathbb{E}(\hat{\boldsymbol \beta}) \neq {\boldsymbol \beta}$ 라는 의미이며 *shrinkage* 는 $\lambda$ 에 의하여 원래 $\mathbb{E}(\boldsymbol \beta)$의 값이 원래 $\beta$ 가 가져야할 값보다 전체적으로 ($L^2$-norm 관점에서!) 작게 추정됨을 의미한다. 요약하면, 능형회귀로 얻은 추정량은 **쓰레기**라는 의미이다. 	

	그렇다면 왜 $\hat{\boldsymbol \beta}=({\bf X}^\top{\bf X}+\lambda {\bf I})^{-1}{\bf X}^\top {\bf y}$ 을 쓰는것일까? $\mathbb{E}(\hat{\boldsymbol \beta})$ 관점에서는 쓰레기가 맞는데 $\mathbb{V}(\hat{\boldsymbol \beta})$ 의 관점에서는 좋은면이 있기 때문이다. 
"""

# ╔═╡ e9e42e74-b295-4cc7-94f4-33b23d4fc751
md"""
> 이전 강의노트에서 *biased* 가 아닌 *unbiased* 라고 되어있었습니다. 제가 설명도 *unbiased* 로 잘못했는데요, *biased*로 이해하고 설명들으시면 됩니다. 혼란을 드려 죄송합니다. 
"""

# ╔═╡ 59e96cb9-f277-4825-a482-5f819905312f
md"""
### B. 추정량의 분산
"""

# ╔═╡ ccd985b7-978d-4b15-a206-f1bf323c38e1
md"""
-- 방법1: 이론적으로 조사하자.
"""

# ╔═╡ d1e558cc-a4d6-469d-878c-79bc69247eb5
md"""
$\mathbb{V}(\hat{\boldsymbol\beta}) =({\bf X}^\top{\bf X}+\lambda {\bf I})^{-1}{\bf X}^\top{\bf X}({\bf X}^\top{\bf X}+\lambda {\bf I})^{-1}\sigma^2$
"""

# ╔═╡ 5d51fb2a-76ba-426b-b6d8-69d8132be55e
inv(X'X+10*I)X'X*inv(X'X+10*I)*σ^2

# ╔═╡ 6592929c-6d2b-4f0d-bbec-229b0d60ef0a
diag(inv(X'X+10*I)X'X*inv(X'X+10*I)'*σ^2)

# ╔═╡ 3e0fb70f-6855-4867-be33-16d20adef786
md"""
*Figure: $\lambda$에 따른 추정량의 분산변화 (이론)*
"""

# ╔═╡ 22a02082-26e7-452d-bfa4-d9ace93fa11b
let 
	Vβ̂(λ) = diag(inv(X'X+λ*I)X'X*inv(X'X+λ*I)*σ^2)
	λs = [1e1,1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11]
	plot(λ -> Vβ̂(λ)[1], λs,xscale= :log10,yscale= :log10, label="V(β̂₁)")
	plot!(λ -> Vβ̂(λ)[2], label="V(β̂₂)")
	plot!(λ -> Vβ̂(λ)[3], label="V(β̂₃)")
end

# ╔═╡ 30e9a8de-b9f1-4d94-8de3-78b3f1e4a475
md"""
-- 방법2: 시뮬레이션으로 알아보자.
"""

# ╔═╡ 95631387-6cd1-4917-b8c9-878cc7c4d837
md"""
*Figure: $\lambda$에 따른 추정량의 분산변화 (시뮬)*
"""

# ╔═╡ ff6c8d44-7fbc-44b5-8c4c-23e76b34ef6d
let 
	function V̂β̂(λ)
		N = 1000
		Y = (β1*X1 + β2*X2 + β3*X3) .+ σ*randn(n,N)
		B̂ = inv(X'X+λ*I)X'Y
		j = ones(N)
		Vβ̂ = (B̂ .- B̂*j/N).^2 *j/N
		return Vβ̂
	end
	λs = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11]
	plot(λ -> V̂β̂(λ)[1],λs,xscale= :log10,yscale= :log10, label="V̂(β̂₁)")
	plot!(λ -> V̂β̂(λ)[2], label="V̂(β̂₂)")
	plot!(λ -> V̂β̂(λ)[3], label="V̂(β̂₃)")
end

# ╔═╡ 8aeae892-a26e-4c4e-bcb4-a10b2c4d710c
md"""
### C. 추정량의 MSE
"""

# ╔═╡ ab512892-c838-4967-a0e3-dfc04e6a9df0
md"""
ref: <https://en.wikipedia.org/wiki/Mean_squared_error#Estimator>
"""

# ╔═╡ 1f92a68d-48cd-46a2-b373-2631cf1b56da
md"""
-- 능형회귀로 얻은 $\hat{\boldsymbol \beta}$ 의 MSE를 비교하여 보자.
"""

# ╔═╡ 70bcd059-954e-4077-8460-ba366a3271ee
md"""
!!! warning "MSE"
	여기에서 MSE는 일반적으로 머신러닝에서 사용하는 $\frac{1}{n}\sum_{i=1}^{n}(y_i-\hat{y}_i)^2$ 의 개념이 아니다. MSE라는 용어를 predictor에 사용할 경우와 estimator에 사용할 경우가 있는데, predictor에 MSE 라는 용어를 사용할 경우에는 $\frac{1}{n}\sum_{i=1}^{n}(y_i-\hat{y}_i)^2$ 를 지칭하지만 estimator에 사용할 경우는 MSE는 아래를 의미한다. 

	$\begin{align}
	\text{MSE}(\hat{\boldsymbol \beta})&=\mathbb{E}_{\boldsymbol \beta}\big[(\hat{\boldsymbol \beta}-{\boldsymbol \beta})^\top(\hat{\boldsymbol \beta}-{\boldsymbol \beta}) \big]\\
	&=\text{tr}\big(\mathbb{V}(\hat{\boldsymbol \beta})\big)+\big(\mathbb{E}(\hat{\boldsymbol \beta})-{\boldsymbol \beta}\big)^\top\big(\mathbb{E}(\hat{\boldsymbol \beta})-{\boldsymbol \beta}\big)\\&=\text{Var}+\text{Bias}^2
	\end{align}$

	일반적으로 ${\boldsymbol \beta}$ 에 대한 추정량이 가졌으면 좋겠는 좋은 성질은 (1) 편향되어있지 않고 (2) 분산이 작은것 인데 MSE는 이 두 가지 기준을 모두 고려한 좋은 평가방법이다. 
"""

# ╔═╡ 2bacdea6-3f94-44ba-bdb3-c549ee36cbd2
md"""
*Figure: $\lambda$에 따른 추정량의 MSE변화 (이론)*
"""

# ╔═╡ d67d0729-a865-402c-9a87-fd5836b7957c
let 
	function MSEβ̂(λ)
		Bias² = (inv(X'X+λ*I)X'X*β-β)'*(inv(X'X+λ*I)X'X*β-β)
		trVar = tr(inv(X'X+λ*I)X'X*inv(X'X+λ*I)*σ^2)
		return Bias² + trVar
	end
	MSEols = (inv(X'X)X'X*β-β)'*(inv(X'X)X'X*β-β) + tr(inv(X'X)*σ^2)
	#MSEols = tr(inv(X'X)*σ^2) # <-- 사실상이거
	@show MSEols
	λs = [1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11]
	plot(λ -> MSEβ̂(λ),λs,xscale= :log10, yscale= :log10, label="MSE -- ridge")
	plot!(λs,fill(MSEols,length(λs)), label="MSE -- ols", linestyle=:dash)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[compat]
CSV = "~0.10.14"
DataFrames = "~1.6.1"
HTTP = "~1.10.6"
Plots = "~1.40.4"
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "54e136119fa801107056dd2c099fbaa1ecd2eaec"

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

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "6c834533dc1fabd820c1db03c839bf97e45a3fab"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.14"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "4b270d6465eb21ae89b732182c20dc165f8bf9f2"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.25.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "b1c55339b7c6c350ee89f2c1604299660525b248"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.15.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "6cbbd4d241d7e6579ab354737f4dd95ca43946e1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.1"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

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
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "9f00e42f8d99fdde64d40c8ea5d14269a2e2c1aa"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.21"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "ddda044ca260ee324c5fc07edb6d7cf3f0b9c350"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.5"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "278e5e0f820178e8a26df3184fcb2280717c79b1"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.5+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

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
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

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

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "e0b5cd21dc1b44ec6e64f351976f961e6f31d6c4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.3"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3da7367955dcc5c54c1ba4d402ccdc09a1a3e046"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+1"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "442e1e7ac27dd5ff8825c3fa62fbd1e86397974b"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.4"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

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

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

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
git-tree-sha1 = "90b4f68892337554d31cdcdbe19e48989f26c7e6"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.3"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

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

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

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

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "5d54d076465da49d6746c647022f3b3674e64156"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.8"
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

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "dd260903fdabea27d9b6021689b3cd5401a57748"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.20.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "52ff2af32e591541550bd753c0da8b9bc92bb9d9"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.7+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

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
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─7ec1510e-130c-11ef-1e0e-518d9cc440ea
# ╟─a366d78e-5826-4a65-a3fe-3469fee80f9a
# ╟─4c8cb126-08bd-478e-8472-fd1547e0d128
# ╟─2d2f8c88-07f9-4e8b-86ae-fc9828aba187
# ╠═128613c3-baba-426a-adbf-8abed179eb49
# ╠═3317fffc-955d-4073-9123-17688a4c6f61
# ╠═18996e8e-c802-4e03-a3c6-c4d83aa5afdc
# ╟─31f2972f-c83f-42ce-999c-842c157899d3
# ╠═b5277d9e-216f-444f-bcae-e7ddf5d2e9b5
# ╠═df19fb9c-8f86-4d55-9899-c868c3512411
# ╠═f9f88913-a23a-4f02-a1cd-f255feb8c834
# ╠═b45868ea-681e-41e3-b9c0-2d0c918e574d
# ╟─d055be6a-aef8-47d0-b0a2-a56b85089bfc
# ╟─8b6bafd7-7dd9-4c2b-a059-d18fc2ad22f8
# ╟─57b42243-1fa1-4a69-b116-044b30d7819b
# ╟─b70e73b7-c2ae-4c94-9f85-ad0195b5858a
# ╟─6601a37b-cb7c-417d-a438-c5c18c19e0cf
# ╟─23d619d0-e790-4ddb-97ec-735e42fa85bc
# ╠═68d0574b-dbd6-4b68-920e-e672fe876ec6
# ╟─ac6df6a5-6742-446a-b119-a0a685de2f3f
# ╠═cc6f7ff7-e7ac-409e-9f63-145c58e4ade5
# ╟─f05be255-fae4-4421-a920-c863d0c389a6
# ╠═cbb62ff1-fe39-468e-91a6-0431192fb25d
# ╟─3bad3692-7f76-494c-93d5-f3951331c74e
# ╟─758668d4-7d98-46db-8f8a-9a81abe6ad25
# ╟─a82b89e1-ac4b-45e0-9e46-db58bc303812
# ╠═a8e9091b-aa4b-49c8-804f-b0db309329a8
# ╟─cbe878f4-07b2-4f21-8a63-7468b2633b6b
# ╟─90ff8925-e0d6-4982-b19d-0f78b5016545
# ╟─85ad9c84-87c8-49d9-9139-63605442aa2d
# ╟─6216c9c1-3d3a-4e2c-9fd7-d18cba3b6147
# ╠═5ba39dcd-1f21-4c83-a420-59de892a8121
# ╟─1d20583d-bddb-4c38-9da5-03d9e4c53f03
# ╠═d3416f40-9387-4a44-986f-8657c8ea1c00
# ╟─324c7de6-c404-470a-afe3-cc756cf9fa94
# ╟─1109294c-c93c-4156-a301-e76260f5b1a4
# ╟─05665ffd-fba5-49c4-bae7-65c6b02ebc64
# ╠═278ef454-7fcc-4155-89f4-6fd93d168feb
# ╟─14f7efbc-448d-44e1-b107-d52f05bc8520
# ╠═c5bdf20e-e89f-46e6-9a9a-2c16a19a54c4
# ╟─21a75194-76de-43bd-abbe-89d3e28e6d8e
# ╠═14e02673-0863-4aa3-9ba6-af0bc909d0cc
# ╟─6d545544-3516-4560-8887-e4413f3c0fc0
# ╟─f5b49064-9f8b-4238-97a4-780b4023a641
# ╟─8a47abe4-22de-48fc-b98c-c4cce53502bb
# ╟─b3d527a8-3151-4f30-8ccd-b53b298855de
# ╠═2fed1588-f9cc-4c6d-a8ff-751352a1fcd0
# ╟─61b295ca-be2d-4326-bff5-0fd860d31919
# ╟─057f12ba-1e54-481a-b9f4-a6e2bcad6805
# ╠═392d159d-c7a4-47f8-bdfd-c28fdfc1e6a6
# ╟─e1fea1fb-7a22-41de-b84a-b02c77c728c9
# ╟─1c21d098-4764-469a-9a93-9c091304b7cd
# ╟─42cd97a0-c604-4656-9d4c-fc4a8f268580
# ╟─cc218bfe-8862-4e48-8302-88b6a17da2a2
# ╟─fbeb2602-7d61-4a6e-9ca6-bf8618721d1d
# ╟─d7206a25-b7f6-4c59-8efc-268f3485fe2d
# ╟─ebb19d9a-0ddc-40a5-9533-79639bba2e82
# ╟─b71607cd-dfb5-4138-9042-fdd1b7b93bba
# ╠═096103d0-b7e3-4f60-b7d9-a69a58b4d14e
# ╟─0401b893-61ed-426a-8f0c-cfa365c56022
# ╠═80b71ff3-3b35-4a6d-bacb-5a69d735dabc
# ╟─86ee8e65-1220-4d8f-8801-a51e1924cd85
# ╠═828c4d10-263b-46d2-b36f-28f6a0f3c612
# ╟─b1266c8b-2be6-4fb0-a791-38d9f2844d0b
# ╟─101cd18a-0b4b-4e1e-8dee-6d687ebedc30
# ╠═c029a13f-ab20-47f0-92e0-edd767302f85
# ╟─0b7a7efb-e3e2-42f5-840c-e8245660e840
# ╟─e9e42e74-b295-4cc7-94f4-33b23d4fc751
# ╟─59e96cb9-f277-4825-a482-5f819905312f
# ╟─ccd985b7-978d-4b15-a206-f1bf323c38e1
# ╟─d1e558cc-a4d6-469d-878c-79bc69247eb5
# ╠═5d51fb2a-76ba-426b-b6d8-69d8132be55e
# ╠═6592929c-6d2b-4f0d-bbec-229b0d60ef0a
# ╟─3e0fb70f-6855-4867-be33-16d20adef786
# ╠═22a02082-26e7-452d-bfa4-d9ace93fa11b
# ╟─30e9a8de-b9f1-4d94-8de3-78b3f1e4a475
# ╟─95631387-6cd1-4917-b8c9-878cc7c4d837
# ╠═ff6c8d44-7fbc-44b5-8c4c-23e76b34ef6d
# ╟─8aeae892-a26e-4c4e-bcb4-a10b2c4d710c
# ╟─ab512892-c838-4967-a0e3-dfc04e6a9df0
# ╟─1f92a68d-48cd-46a2-b373-2631cf1b56da
# ╟─70bcd059-954e-4077-8460-ba366a3271ee
# ╟─2bacdea6-3f94-44ba-bdb3-c549ee36cbd2
# ╠═d67d0729-a865-402c-9a87-fd5836b7957c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
