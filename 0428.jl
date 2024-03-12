### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ a6609191-8981-4d00-a5ef-b311d8791a55
using PlutoUI,Plots,Distributions,Random,ForwardDiff

# ╔═╡ 771b13aa-2bdd-400a-87e3-cc5bae383cfe
md"""
# 4월28일: 통계전산 2022 중간고사 풀이
"""

# ╔═╡ d1feebb2-87f1-4cdb-a5ec-9d4755673682
html"""
<div style="display: flex; justify-content: center;">
<div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
<iframe src=
"
https://www.youtube.com/embed/playlist?list=PLQqh36zP38-wVkXbNqaC9MDXhrRIKVfXu
"
width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ a8260ce7-b754-4423-a626-a3db25cf12ed
Plots.plotly()

# ╔═╡ 0ccda176-3037-4c56-8d22-432597c8014b
md"""
## 1. 줄리아문법 (5점)
"""

# ╔═╡ 6f27a359-bf9e-48b8-a4fd-26c9d499ba8e
md"""
(1) ForwardDiff 를 이용하여 $f(x)=\frac{1}{\sqrt{2\pi}}e^{-\frac{1}{2}x^2}$ 일때 $f'(0)$ 의 값을 계산하라. 
"""

# ╔═╡ af42267f-4682-48a9-af8e-8c1419bfbfae
let 
	f(x) = 1/(√(2π)) * exp(-0.5*x^2)
	ForwardDiff.derivative(f,0)
end

# ╔═╡ 5f7c6cfc-878c-4266-953b-21419ba07686
md"""
(2) randn(100) 을 이용하여 Array를 생성하고 양수를 counting 하라. 
"""

# ╔═╡ 7200e56a-e6b5-42c7-bc8d-574ef4cd6b55
sum(randn(100) .> 0)

# ╔═╡ 92dfb9c1-f1e9-4de7-8e91-d59384ebfb86
md"""
(3) 함수 $f(x)=\frac{e^x}{1+e^x}$ 를 (-5,5)의 범위에서 plot 하라. 
"""

# ╔═╡ 5342b97d-b7a6-4552-9141-6b16155c01f3
let 
	f(x)= exp(x) / (1+exp(x))
	plot(f,-5,5)
end 

# ╔═╡ 069d6e5e-f5a1-47c9-b4f3-c9de6233cff0
md"""
(4) 함수 $f(x) = \max(0,x)$ 를 (-5,5)의 범위에서 plot 하라. 
"""

# ╔═╡ c2d413c7-d19d-4868-aa35-6b1d86d85741
let 
	f(x)=max(0,x)
	plot(f,-5,5) 
end

# ╔═╡ ab443d25-6865-48e5-a5ff-998d8c743890
md"""
(5) 아래와 같은 수열을 100항까지 생성하라. 
- 1, 1+2, 1+2+3,... 
"""

# ╔═╡ 96a7319f-2b02-4498-91c9-0da225940e53
[sum(1:i) for i in 1:100]

# ╔═╡ 018f278c-c373-11ec-1e5b-57380df5b3ff
md"""
## 2. 몬테카를로 적분 (10점)
몬테카를로 방법을 활용하여 아래의 적분을 계산하라. 
"""

# ╔═╡ ea1bfd46-418e-415f-8f6b-8da07287ab82
md"""
(1) $\int_{0}^{1}\sin(x)dx$ 
"""

# ╔═╡ 68e72126-ba72-4afd-adac-53e80b289fa3
let 
	X = rand(1000)
	mean(sin.(X))
end 

# ╔═╡ dc51c540-d936-48ae-9144-20ce3717ae9f
md"""
(2) $\int_{0}^{\infty} \big(\frac{1}{2}\big)^xdx$
"""

# ╔═╡ cedec46d-b5fa-42ea-b45d-06c0b8707f3b
let 
	X = rand(Exponential(1),1000) 
	mean(@. (0.5)^X * exp(X)) 
end 

# ╔═╡ 333fc0c4-648a-4c6d-9803-d5e0b5eb6640
md"""
- note: $\int_0^{\infty}(1/2)^xdx = \int_0^{\infty}(1/2)^xe^x e^{-x}dx=E\big[(0.5)^X e^X\big], \quad X \sim Exp(1)$
"""

# ╔═╡ dd913eda-4442-4948-83e6-e8b1da7d3646
md"""
(3) $\int_{3}^{\infty}e^{-x}dx$
"""

# ╔═╡ 6ee56031-39ac-4f63-8d2b-e4b218338d33
let 
	X = rand(Exponential(1),1000)
	mean(X .>3)
end 

# ╔═╡ c798d519-3afb-4265-b69f-a884cfcb7a55
md"""
(4) $\int_{10000}^{\infty}e^{-x}dx / \int_{9999}^{\infty}e^{-x}dx$
"""

# ╔═╡ ab5287aa-3b85-47e1-b271-f861854c170a
md"""
- note: $\int_{10000}^{\infty}e^{-x}dx / \int_{9999}^{\infty}e^{-x}dx=\int_{1}^{\infty}e^{-x}dx / \int_{0}^{\infty}e^{-x}dx=\int_{1}^{\infty}e^{-x}dx$
"""

# ╔═╡ 4bb3343f-a789-4a5b-9b11-77805fce7f91
let 
	X= rand(Exponential(1),1000)
	mean(X .> 1)
end

# ╔═╡ 782a4764-3709-4e86-b990-acc3e10ae2e3
md"""
(5) $\int_{-1.96}^{1.96} \frac{1}{\sqrt{2\pi}}e^{-x^2/2}dx$ 
"""

# ╔═╡ d29923cf-298c-49d2-9bc2-ffaa91c54f8b
let 
	X =rand(Normal(0,1),1000)
	mean( -1.96 .< X .< 1.96)
end 

# ╔═╡ aa306e3f-6e9a-413c-8e75-22f49fbc150a
md"""
## 3. 난수생성 (30점)
"""

# ╔═╡ b80b9299-7647-4296-b70a-a39e47d75d16
md"""
(1) $U(0,1)$ 로 부터 $X_1,\dots, X_{1000} \overset{iid}{\sim} bin(5,0.4)$를 생성하라. 
"""

# ╔═╡ aa9a809f-842f-4dc2-aa65-6084a8cf380a
[sum(rand(5) .< 0.4) for i in 1:1000] 

# ╔═╡ b2d9e3bf-25c0-498b-8439-41c0a4f8f615
md"""
(2) $Poi(1)$ 로 부터 $X_1,\dots, X_{1000} \overset{iid}{\sim} Poi(10)$을 생성하라. 
"""

# ╔═╡ ce3174aa-2a20-4893-86ec-78bbd5e3a396
[sum(rand(Poisson(1),10)) for i in 1:1000]

# ╔═╡ 15926539-ccac-4796-bafb-4d828145037e
md"""
(3) $Ber(0.01)$ 로 부터 $X_1, \dots, X_{1000} \overset{iid}{\sim} Poi(7)$을 근사적으로 생성하라. 
"""

# ╔═╡ 62f75765-d05c-49bf-9085-d54d56daedc4
[sum(rand(Bernoulli(0.01),700)) for i in 1:1000]

# ╔═╡ 41239d55-e7d7-4727-b409-3fbc01fc2fc5
md"""
(4) $\Gamma(\frac{1}{2},4)$ 로 부터 $X_1,\dots,X_{1000} \overset{iid}{\sim} \chi^2(56)$을 생성하라.
"""

# ╔═╡ 483db57f-d04b-4edb-8921-712d47ac414c
[sum(rand(Gamma(0.5,4),56)./2) for i in 1:1000]

# ╔═╡ 627b5bcf-d3a4-412e-aee4-a9a5df0bb1b0
md"""
(5) $\Gamma(\frac{1}{2},4)$ 로 부터 $X_1,\dots, X_{1000} \overset{iid}{\sim} Exp(2)$를 생성하라. 
"""

# ╔═╡ 2792df06-4618-4d08-a797-b8f57455b6bd
[rand(Gamma(0.5,4))/2+rand(Gamma(0.5,4))/2 for i in 1:1000] # 평균이2인 지수분포

# ╔═╡ 63cde8af-626c-4814-9e93-f4850178d3b2
[(rand(Gamma(0.5,4))/2+rand(Gamma(0.5,4))/2)/4 for i in 1:1000] # 평균이1/2인 지수분포

# ╔═╡ 253c7b93-afa6-49e8-b615-1aa64a22b168
md"""
(6) $Exp(2)$ 로 부터 $X_1,\dots,X_{1000} \overset{iid}{\sim} \Gamma(10,4)$를 생성하라.
"""

# ╔═╡ bbbca719-599b-4bf4-b9de-40a7565421dd
[sum(rand(Exponential(2),10) .*2)  for i in 1:1000] # 평균이2인 지수분포 

# ╔═╡ e4f4d6e4-847d-4295-bf9c-cdfc261657f3
[sum(rand(Exponential(0.5),10) .*8)  for i in 1:1000] # 평균이1/2인 지수분포 

# ╔═╡ 501ce3c3-891c-48d6-a9a5-c035072719fb
md"""
(7) $\chi^2(1)$ 로 부터 $X_1,\dots,X_{1000} \overset{iid}{\sim}Exp(10)$을 생성하라. 
- 여기에서 $Exp(10)$는 평균이 1/10인 지수분포이다. 즉 $X\sim Exp(10)$이면 $X$의 pdf는 $f(x)=\frac{1}{10}e^{-x/10}$ 으로 주어진다. 
"""

# ╔═╡ f5aa6294-de7b-4821-acd2-1ea61745cbe8
[(rand(Chisq(1)) + rand(Chisq(1)))*5 for i in 1:1000]  # 평균이 10인 지수분포 

# ╔═╡ 79a209e7-d885-46f7-8b8d-b530c0cb8150
[(rand(Chisq(1)) + rand(Chisq(1)))/20 for i in 1:1000]  # 평균이 1/10인 지수분포 

# ╔═╡ 908f12ab-790c-4a64-b701-0cb6021c4fc4
md"""
(8) $U(0,1)$ 로 부터 $X_1, \dots, X_{1000} \overset{iid}{\sim} \Gamma(10,2)$을 생성하라.
"""

# ╔═╡ 5985edc7-3718-44f8-b13b-62614d04d5b1
[sum(-log.(1 .- rand(10)))*2 for i in 1:10000]

# ╔═╡ 3353dd22-4344-4f09-80ef-3c244d4254bf
md"""
(9) $U(0,1)$ 로 부터 $X_1,\dots,X_{1000} \overset{iid}{\sim} N(0,1)$ 을 생성하라. 
"""

# ╔═╡ 11b0c663-d04b-4de8-a0cf-4272ea9e8ee1
let 
	U₁ = rand(1000)
	U₂ = rand(1000)
	@. √(-2log(1-U₁))*cos(2π*U₂)
end

# ╔═╡ e8ded97f-fa69-4f5f-a85c-74faa4911eba
md"""
(10) $N(0,1)$ 로 부터 $X_1,\dots,X_{1000} \overset{iid}{\sim} Exp(1)$ 을 생성하라. 
"""

# ╔═╡ bd0d44d6-21bd-4346-a1ec-29bedb81237a
[sum(rand(Normal(0,1),2).^2)/2 for i in 1:1000]

# ╔═╡ 93dfe906-ceb8-477e-9146-dc952e2052c0
md"""
## 4. Delta Method (30점)
"""

# ╔═╡ aafae26c-8e6d-483a-b8ff-c18696c61744
md"""
 $X_1,\dots,X_n \sim Ber(p)$이라고 하자. 모수의 추정치를 $\hat{p}=\bar{X}_n=\frac{1}{n}\sum_{i=1}^{n}X_n$라고 하자. 아래의 물음에 답하라. 

(1) 평균은 $\hat{p}$로 분산은 $\hat{p}(1-\hat{p})$로 추정한다고 하자. $p=0.4$일 경우 $\mu$와 $\sigma^2$에 대한 95% 점근적 신뢰구간을 시뮬레이션과 delta method를 이용하여 구하고 비교하라. 

"""

# ╔═╡ 144356a9-b763-4365-9c89-6e295d2c650d
# p̂ 신뢰구간 시뮬 
let 
	n=1000
	p=0.4 
	estimators = [mean(rand(Bernoulli(p),n)) for i in 1:10000]
	quantile(estimators,0.025), quantile(estimators,0.975)
end

# ╔═╡ 94891d36-63d5-4e67-93d2-44d3d675b408
# p̂ 신뢰구간 이론 
let 
	n=1000
	p=0.4 
	adist = Normal(p,√(p*(1-p)/n)) ## CLT 
	quantile(adist,0.025), quantile(adist,0.975)
end 

# ╔═╡ 9ed1cd40-a1db-49a8-b10c-f346cfb816d0
# p̂(1-p̂) 신뢰구간 시뮬 
let 
	n=1000
	p=0.4 
	g(x)=x*(1-x)
	estimators = [mean(rand(Bernoulli(p),n)) .|> g for i in 1:10000]
	quantile(estimators,0.025), quantile(estimators,0.975)
end

# ╔═╡ bd7ee6f2-4be0-4fde-a303-822a9095d88f
# p̂(1-p̂) 신뢰구간 이론 1
let 
	n=1000
	p=0.4 
	g(x)=x*(1-x)
	adist = Normal(g(p),√(p*(1-p)*(1-2p)^2/n)) ## CLT + delta method
	quantile(adist,0.025), quantile(adist,0.975)
end 

# ╔═╡ 5eed3134-8b1d-4975-97a4-2b6d8a0a99bc
# p̂(1-p̂) 신뢰구간 이론 2
let 
	n=1000
	p=0.4 
	g(x)=x*(1-x)
	g′ = x-> ForwardDiff.derivative(g,x)
	adist = Normal(g(p),√(p*(1-p)/n)*g′(p)) ## CLT + delta method
	quantile(adist,0.025), quantile(adist,0.975)
end 

# ╔═╡ 045e8fde-29f2-4cf4-97f7-c53af4c7be61
md"""
(2) $g(x)=\sin^{-1}(\sqrt{x})$ 일때 $p=0.4$에 대한 $g(\hat{p})$의 95% 점근적 신뢰구간을 시뮬레이션과 delta method를 이용하여 구하고 비교하라. 
- 힌트: $\frac{d}{du}\sin^{-1}(u)=1/\sqrt{1-u^2}$.
"""

# ╔═╡ 9c98fb88-99cd-4ffc-9e6f-0656734a501e
# g(p̂) 신뢰구간 시뮬 
let 
	n=1000
	p=0.4 
	g(x)=asin(√x)
	estimators = [mean(rand(Bernoulli(p),n)) .|> g for i in 1:10000]
	quantile(estimators,0.025), quantile(estimators,0.975)
end

# ╔═╡ 7f4e326f-624e-4e2d-a814-b9ad5445ce63
# p̂(1-p̂) 신뢰구간 이론 2
let 
	n=1000
	p=0.4 
	g(x)=asin(√x)
	g′ = x-> ForwardDiff.derivative(g,x)
	adist = Normal(g(p),√(p*(1-p)/n)*g′(p)) ## CLT + delta method
	quantile(adist,0.025), quantile(adist,0.975)
end 

# ╔═╡ fcdb21a6-ebbc-4797-9239-87cdc0eeaf73
md"""
(3) 임의의 $p$에 대한 $g(\hat{p})$에 대한 분산을 시뮬레이션과 delta method를 이용하여 구하고 비교하라.  $p$에 따른 분산변화가 있는가? 
"""

# ╔═╡ fc0ec748-0a82-4cb1-98e3-0b709a320ba3
# g(p̂) 신뢰구간 시뮬 
let 
	n=1000
	p=0.897
	g(x)=asin(√x)
	estimators = [mean(rand(Bernoulli(p),n)) .|> g for i in 1:10000]
	var(estimators)
end

# ╔═╡ 1ee351a8-bc0b-4ace-bb6b-9d2f5e5e32b7
# p̂(1-p̂) 신뢰구간 이론 2
let 
	n=1000
	p=0.345456
	g(x)=asin(√x)
	g′ = x-> ForwardDiff.derivative(g,x)
	σ=√(p*(1-p)/n)*g′(p)
	σ^2 # 1/(4n)
end 

# ╔═╡ 9a2ccacf-c2dc-49a2-93ad-7b52f2540c08
md"""
`-` $p$ 에 따른 분산변화가 없다.
"""

# ╔═╡ 3cbaa153-6a09-414f-beaf-40d50faa01b4
md"""
## 5. 가설검정과 신뢰구간 (25점)

(1) $X_1\dots X_n \overset{iid}{\sim} N(\mu,1)$ 에서 아래를 test하고 싶다고 하자. 

$H_0: \mu=0$

$H_1:\mu>0$

이를 위하여 검정통계량으로 $\bar{x}=\frac{1}{30}\sum_{i=1}^{30}x_i=0.01$을 얻었다고 하자. $p$-value의 이론값과 시뮬레이션 값을 계산하고 비교하여라.  
"""

# ╔═╡ 7012b1fc-ae96-40e8-a1f0-848c34fb9174
## 시뮬 
([mean(rand(Normal(0,1),30)) for i in 1:10000] .> 0.01) |> mean

# ╔═╡ 3f31cc87-7962-49c8-9b77-9a0d44e4cd09
## 이론
1-cdf(Normal(0,1/√30),0.01)

# ╔═╡ 20a6428d-de79-4c2b-bfeb-823ba56b94cf
md"""
(2) $X_1,\dots,X_n \overset{iid}{\sim} N(\mu,15^2)$ 이라고 하자. 구간 $(\bar{x}-1,\bar{x}+1)$이 $\mu$에 대한 95% 신뢰구간으로 주어지기 위한 표본크기 $n$을 구하여라. (95%이상의 신뢰구간을 만드는 가장 작은 $n$을 구할 것) 
"""

# ╔═╡ 119cfb0c-98f3-4ae6-a8b8-e29bd4d85e9b
let 
	n=865
	quantile(Normal(0,15/√n),0.025), quantile(Normal(0,15/√n),0.975)
end 

# ╔═╡ 24d06f21-3ccc-4be7-ad6d-a13aa39b5f59
md"""
--- 
**(3)-(5)**

 $X_1, \dots, X_n \sim N(\mu_x,1)$ 이고 $Y_1, \dots, Y_n \sim N(\mu_y,1)$이라고 하자. 아래를 test 하고 싶다고 하자. 

$H_0: (\mu_x,\mu_y)=(0,0)$

$H_1: (\mu_x,\mu_y)\neq (0,0)$

이를 위한 검정통계량으로 $r^2=\bar{x}^2 + \bar{y}^2$ 을 사용한다고 하자. 3개의 샘플을 관찰하여 

-  $\bar{x}=\frac{1}{3}\sum_{i=1}^{3}x_i=0.8$ 
-  $\bar{y}=\frac{1}{3}\sum_{i=1}^{3}y_i=0.7$ 

을 얻었다고 하자.
"""

# ╔═╡ 3fec2a40-315e-40dc-b77b-6f3ae1b8aaa8
md"""
(3) 90%의 유의수준에서 귀무가설을 기각 (= $H_1$을 채택) 하기 위한 가장 작은 $r$의 이론값과 시뮬레이션 값을 각각 구하고 비교하라. 
"""

# ╔═╡ a322512e-7dba-4481-b1cb-c60fcbfc9bec
# 시뮬
let
	dst = Normal(0,1/√3)
	r= .√(rand(dst,1000).^2 + rand(dst,1000).^2)
	quantile(r,0.9)
end

# ╔═╡ 064a47d3-7753-42f9-b919-dcb4333858f9
# 이론
let
	quantile(Chisq(2),0.9) # (√3X̄)^2+(√3Ȳ)^2 = 3R² ~ Χ²(2)
	√(quantile(Chisq(2),0.9)/3)
end

# ╔═╡ 56dbab84-ca18-4fcf-9ea4-86403215a7c3
md"""
(4) 주어진 검정통계량의 $p$-value의 이론값과 시뮬레이션값을 각각 구하고 비교하라. 
"""

# ╔═╡ db58ac27-9b24-40da-bbea-5a5ee0a0be30
# 시뮬
let 
	dst = Normal(0,1/√3)
	r= .√(rand(dst,1000).^2 + rand(dst,1000).^2)
	(r.^2 .> 0.8^2+0.7^2) |> mean
end

# ╔═╡ 08fa30d5-70b9-46e0-9515-f5d32e85b3a1
# 이론 
let 
	1-cdf(Chisq(2), 3(0.8^2+0.7^2)) # (√3X̄)^2+(√3Ȳ)^2 = 3R² ~ Χ²(2)
end

# ╔═╡ e50b8476-be2f-40fb-b054-2866fbefc3db
md"""
(5) 샘플수가 30개 일 경우 $\bar{x}=\frac{1}{30}\sum_{i=1}^{30}x_i=0.8$ 과 $\bar{y}=\frac{1}{30}\sum_{i=1}^{30}y_i=0.7$ 을 얻었다면 $p$-value는 어떻게 되는가? (4)의 $p$-value와 비교하여라. 
-  $n=30$일 경우 $p$-value는 이론, 시뮬레이션중 어떠한 하나의 방식으로 계산하여도 무방함. 
"""

# ╔═╡ d51b8ca8-d2b7-4e70-843d-aab78127f2e4
# 시뮬
let 
	dst = Normal(0,1/√30)
	r= .√(rand(dst,1000).^2 + rand(dst,1000).^2)
	(r.^2 .> 0.8^2+0.7^2) |> mean
end

# ╔═╡ 0e28a45a-602f-4c94-85ec-261680507863
# 이론 
let 
	1-cdf(Chisq(2), 30(0.8^2+0.7^2)) # (√30X̄)^2+(√30Ȳ)^2 = 30R² ~ Χ²(2)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
Distributions = "~0.25.53"
ForwardDiff = "~0.10.27"
Plots = "~1.27.6"
PlutoUI = "~0.7.38"
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
git-tree-sha1 = "9950387274246d08af38f6eef8cb5480862a435f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.14.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "12fc73e5e0af68ad3137b886e3f7c1eacfca2640"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.17.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

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
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

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

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "dd933c4ef7b4c270aacd4eb88fa64c147492acf0"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.10.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "5a4168170ede913a2cd679e53c2123cb4b889795"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.53"

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

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

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

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "34e6147e7686a101c245f12dba43b743c7afda96"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.27"

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
git-tree-sha1 = "af237c08bda486b74318c8070adb96efa6952530"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.2"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "cd6efcf9dc746b06709df14e462f0a3fe0786b1e"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.2+0"

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
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

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
git-tree-sha1 = "91b5dcf362c5add98049e6c29ee756910b03051d"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.3"

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
git-tree-sha1 = "6f14549f7760d84b2db7a9b10b88cd3cc3025730"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.14"

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
git-tree-sha1 = "a970d55c2ad8084ca317a4658ba6ce99b7523571"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.12"

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

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "e8185b83b9fc56eb6456200e873ce598ebc7f262"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.7"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "3b429f37de37f1fc603cc1de4a799dc7fbe4c0b6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.0"

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
git-tree-sha1 = "6f2dd1cf7a4bbf4f305a0d8750e351cb46dfbe80"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.27.6"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "670e559e5c8e191ded66fa9ea89c97f10376bb4c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.38"

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

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

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

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

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
git-tree-sha1 = "5ba658aeecaaf96923dce0da9e703bd1fe7666f9"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.4"

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

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5950925ff997ed6fb3e985dcce8eb1ba42a0bbe7"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.18"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "57617b34fa34f91d536eb265df67c2d4519b8b98"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.5"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

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

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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
# ╟─771b13aa-2bdd-400a-87e3-cc5bae383cfe
# ╟─d1feebb2-87f1-4cdb-a5ec-9d4755673682
# ╠═a6609191-8981-4d00-a5ef-b311d8791a55
# ╠═a8260ce7-b754-4423-a626-a3db25cf12ed
# ╟─0ccda176-3037-4c56-8d22-432597c8014b
# ╟─6f27a359-bf9e-48b8-a4fd-26c9d499ba8e
# ╠═af42267f-4682-48a9-af8e-8c1419bfbfae
# ╟─5f7c6cfc-878c-4266-953b-21419ba07686
# ╠═7200e56a-e6b5-42c7-bc8d-574ef4cd6b55
# ╟─92dfb9c1-f1e9-4de7-8e91-d59384ebfb86
# ╠═5342b97d-b7a6-4552-9141-6b16155c01f3
# ╟─069d6e5e-f5a1-47c9-b4f3-c9de6233cff0
# ╠═c2d413c7-d19d-4868-aa35-6b1d86d85741
# ╟─ab443d25-6865-48e5-a5ff-998d8c743890
# ╠═96a7319f-2b02-4498-91c9-0da225940e53
# ╟─018f278c-c373-11ec-1e5b-57380df5b3ff
# ╟─ea1bfd46-418e-415f-8f6b-8da07287ab82
# ╠═68e72126-ba72-4afd-adac-53e80b289fa3
# ╟─dc51c540-d936-48ae-9144-20ce3717ae9f
# ╠═cedec46d-b5fa-42ea-b45d-06c0b8707f3b
# ╟─333fc0c4-648a-4c6d-9803-d5e0b5eb6640
# ╟─dd913eda-4442-4948-83e6-e8b1da7d3646
# ╠═6ee56031-39ac-4f63-8d2b-e4b218338d33
# ╟─c798d519-3afb-4265-b69f-a884cfcb7a55
# ╟─ab5287aa-3b85-47e1-b271-f861854c170a
# ╠═4bb3343f-a789-4a5b-9b11-77805fce7f91
# ╟─782a4764-3709-4e86-b990-acc3e10ae2e3
# ╠═d29923cf-298c-49d2-9bc2-ffaa91c54f8b
# ╟─aa306e3f-6e9a-413c-8e75-22f49fbc150a
# ╟─b80b9299-7647-4296-b70a-a39e47d75d16
# ╠═aa9a809f-842f-4dc2-aa65-6084a8cf380a
# ╟─b2d9e3bf-25c0-498b-8439-41c0a4f8f615
# ╠═ce3174aa-2a20-4893-86ec-78bbd5e3a396
# ╟─15926539-ccac-4796-bafb-4d828145037e
# ╠═62f75765-d05c-49bf-9085-d54d56daedc4
# ╟─41239d55-e7d7-4727-b409-3fbc01fc2fc5
# ╠═483db57f-d04b-4edb-8921-712d47ac414c
# ╟─627b5bcf-d3a4-412e-aee4-a9a5df0bb1b0
# ╠═2792df06-4618-4d08-a797-b8f57455b6bd
# ╠═63cde8af-626c-4814-9e93-f4850178d3b2
# ╟─253c7b93-afa6-49e8-b615-1aa64a22b168
# ╠═bbbca719-599b-4bf4-b9de-40a7565421dd
# ╠═e4f4d6e4-847d-4295-bf9c-cdfc261657f3
# ╟─501ce3c3-891c-48d6-a9a5-c035072719fb
# ╠═f5aa6294-de7b-4821-acd2-1ea61745cbe8
# ╠═79a209e7-d885-46f7-8b8d-b530c0cb8150
# ╟─908f12ab-790c-4a64-b701-0cb6021c4fc4
# ╠═5985edc7-3718-44f8-b13b-62614d04d5b1
# ╟─3353dd22-4344-4f09-80ef-3c244d4254bf
# ╠═11b0c663-d04b-4de8-a0cf-4272ea9e8ee1
# ╟─e8ded97f-fa69-4f5f-a85c-74faa4911eba
# ╠═bd0d44d6-21bd-4346-a1ec-29bedb81237a
# ╟─93dfe906-ceb8-477e-9146-dc952e2052c0
# ╟─aafae26c-8e6d-483a-b8ff-c18696c61744
# ╠═144356a9-b763-4365-9c89-6e295d2c650d
# ╠═94891d36-63d5-4e67-93d2-44d3d675b408
# ╠═9ed1cd40-a1db-49a8-b10c-f346cfb816d0
# ╠═bd7ee6f2-4be0-4fde-a303-822a9095d88f
# ╠═5eed3134-8b1d-4975-97a4-2b6d8a0a99bc
# ╟─045e8fde-29f2-4cf4-97f7-c53af4c7be61
# ╠═9c98fb88-99cd-4ffc-9e6f-0656734a501e
# ╠═7f4e326f-624e-4e2d-a814-b9ad5445ce63
# ╟─fcdb21a6-ebbc-4797-9239-87cdc0eeaf73
# ╠═fc0ec748-0a82-4cb1-98e3-0b709a320ba3
# ╠═1ee351a8-bc0b-4ace-bb6b-9d2f5e5e32b7
# ╟─9a2ccacf-c2dc-49a2-93ad-7b52f2540c08
# ╟─3cbaa153-6a09-414f-beaf-40d50faa01b4
# ╠═7012b1fc-ae96-40e8-a1f0-848c34fb9174
# ╠═3f31cc87-7962-49c8-9b77-9a0d44e4cd09
# ╟─20a6428d-de79-4c2b-bfeb-823ba56b94cf
# ╠═119cfb0c-98f3-4ae6-a8b8-e29bd4d85e9b
# ╟─24d06f21-3ccc-4be7-ad6d-a13aa39b5f59
# ╟─3fec2a40-315e-40dc-b77b-6f3ae1b8aaa8
# ╠═a322512e-7dba-4481-b1cb-c60fcbfc9bec
# ╠═064a47d3-7753-42f9-b919-dcb4333858f9
# ╟─56dbab84-ca18-4fcf-9ea4-86403215a7c3
# ╠═db58ac27-9b24-40da-bbea-5a5ee0a0be30
# ╠═08fa30d5-70b9-46e0-9515-f5d32e85b3a1
# ╟─e50b8476-be2f-40fb-b054-2866fbefc3db
# ╠═d51b8ca8-d2b7-4e70-843d-aab78127f2e4
# ╠═0e28a45a-602f-4c94-85ec-261680507863
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
