### A Pluto.jl notebook ###
# v0.18.2

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

# ╔═╡ 9d2a3140-7c5e-4439-93e0-1e02e8417daf
using PlutoUI,Distributions,Plots

# ╔═╡ b8113ace-bf87-11ec-332c-d57eb9d49518
md"""
# 4월19일 
"""

# ╔═╡ 1cfafdd4-86ff-4e2f-aa50-e9e78f4f4090
html"""
<div style="display: flex; justify-content: center;">
<div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
<iframe src=
"
https://www.youtube.com/embed/playlist?list=PLQqh36zP38-wxN8yuXPs7sk0zv1VvNJyB
"
width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ e72224fa-ebba-48e4-b702-b2d0329610a6
md"""
## usings 
"""

# ╔═╡ 7fd8042e-9af6-4f6c-b25c-ab8ed5f0af93
PlutoUI.TableOfContents()

# ╔═╡ aa3f4785-8e3f-4841-b3cd-fb1077792992
Plots.plotly()

# ╔═╡ fe50793a-c255-4cc6-81fe-21d729262786
md"""
## 카이제곱분포: $X\sim \chi^2(k)$
"""

# ╔═╡ 823932ac-0546-427f-8d1a-bb1a452ef29c
md"""
### motive 
"""

# ╔═╡ 019f20b3-7920-40e9-b470-4f6392b1c4ac
md"""
(예제) $X_i \overset{iid}{\sim} N(7,\sigma^2)$일때 아래를 test하고 싶다고 하자. 

$H_0: \sigma^2 =4$ 
$H_1: \sigma^2 <4$

30개의 샘플을 확보하여 $\xi=\frac{1}{30}\sum_{i=1}^{30}(x_i-7)^2$를 계산하였으며 계산결과 $\xi=2.72$가 나왔다고 하자. $p$-value를 구하여라. 
"""

# ╔═╡ 930005f3-d5dc-4642-ab96-04f9c74a1c8f
md"""
(풀이1)
"""

# ╔═╡ 0fd62af9-1519-4a64-a860-99e862087320
ξs = [rand(Normal(7,2),30) |> x-> sum((x.-7).^2)/30 for i in 1:14000605]

# ╔═╡ 410876fb-e8d3-413b-b278-9c0d7e76b9a4
mean(ξs.< 2.72)

# ╔═╡ 01e17c60-fd70-4767-b51b-cfd7d542df9f
md"""
- 생각보다 나올법한 확률이다. 
"""

# ╔═╡ 90306510-4d75-4189-89f1-97d3a1f83f21
md"""
`-` 분산이 4가 아닌것 같다라고 주장하기 위해서는 (95%)
"""

# ╔═╡ 239b3846-48df-40b1-b223-fbffa87cd5a4
quantile(ξs,0.05)

# ╔═╡ 076614a8-3c5a-47cc-bcc8-c4074e9776fe
md"""
(풀이2) $\xi=\frac{1}{30}\sum_{i=1}^{30}(x_i-7)^2$는 어떠한 분포A에서 발생한 샘플이라고 볼 수 있다. 그 A의 분포를 이론적으로 잡아보자. 
- 이론: $\sum_{i=1}^{n}Z_i^2\sim \chi^2(n)$, where $Z_i\overset{iid}{\sim} N(0,1)$.
- 관찰: 우리예제의 경우에는 $H_0$가 참이라는 가정하에 $\sum_{i=1}^{30}(\frac{X_i-7}{2})^2\sim \chi^2(30)$
- 주장: $\xi\times 30/4$는 $\chi^2(30)$에서 뽑힌 샘플이다. 
"""

# ╔═╡ 528f54ca-8c1a-4e99-84bb-e1c91004057f
ξ = 2.72 

# ╔═╡ f7274989-9cd5-48dd-b501-91719cc07703
ξ*30/4 

# ╔═╡ 170bf0c7-3055-4c7c-95e2-e2d3d9c3244c
cdf(Chisq(30),ξ*30/4) # 시뮬레이션값과 비슷 

# ╔═╡ 69d57189-a1d8-4574-9714-d5524ff94c55
md"""
`-` 분산이 4는 아닌것같다라고 주장하려면? (95%)
"""

# ╔═╡ 1912ea52-b536-4200-91a1-d9580e635b93
quantile(Chisq(30),0.05)

# ╔═╡ c8a2897c-1381-4828-8d87-d9e3140788e6
md"""
-  $\xi \times 30/4 = 18.49266098195347$ 라는 말이니까 
"""

# ╔═╡ 1d5ced10-a8b3-4e2c-9f0f-e2b7d274c8c3
quantile(Chisq(30),0.05) / 30 *4

# ╔═╡ 3dd308dd-3f64-4a1c-a0e5-96d88d1a68c7
md"""
### 참고: 검정의 형식 논리 
"""

# ╔═╡ abf6b2ce-b5b0-4921-8edd-bda6b38ff364
md"""
`-` 검정을 진행하는 방법은 아래와 같다. 
- 기: 누군가가 (혹은 세상이) $H_0$가 참이라고 주장한다. 나는 $H_1$이 참인것 같다. 
- 승: 누군가와 (혹은 세상과) 싸우기 위하여 샘플을 수집하고 검정통계량을 구한다. 
- 전: 검정통계량의 분포를 잡아내서 $p$-value를 계산한다. 이 $p$-value는 "니가 틀렸겠지"라는 주장에 대한 카운터. 
- 결: $H_0$가 참일지 $H_1$이 참일지 판단. 절대적인 기준은 없음. (하지만 굉장히 보수적인 사람이라도 $p$-value가 0.05보다 작으면 $H_1$이 참이라고 인정)
"""

# ╔═╡ 67d81f4b-4b16-4b0b-aacb-1d38a6e1dbc5
md"""
`-` 포인트는 검정통계량의 분포를 잡아내는 것인데, $H_0$가 참이라는 전제하에 시뮬레이팅 해도 되고 이론적인 분포를 손으로 유도해도 된다. 
- 당연히 컴퓨터가 없던 시절에는 시뮬레이팅이 불가능했으므로 "이론적으로유도+통계표(?)"를 이용해서 $p$-value를 계산해야 했다. 
"""

# ╔═╡ 9061c67f-5039-4ee2-a749-d51aa8034320
md"""
`-` 다양한 분포를 공부하는 이유? 검정통계량의 이론적 분포를 잡아내기 위해서! + α
- 카이제곱분포를 왜 공부해야? 정규분포를 따르는 샘플의 분산을 test하기 위해서! + α 
"""

# ╔═╡ 7fa5ef2f-6aa2-4c00-bf4b-fc47c409dd10
md"""
### 카이제곱분포 요약 
- X의의미: 서로 독립인 표준정규분포의 제곱합 
- X의범위: $x\in (0,\infty)$
- 파라메터의 의미: $k$는 자유도, 표준정규분포 제곱을 몇개 합쳤는지.. 
- 파라메터의 범위: $k=1,2,3,4,\dots$ 
- pdf: 
- mgf: 
- E(X): $k$
- V(X): $2K$
"""

# ╔═╡ a46d857a-6328-43c2-8aff-7155462c2119
md"""
### 대의적정의 
"""

# ╔═╡ 27a1f077-e674-40ab-889c-06663a133c8a
md"""
`-` $X \sim \chi^2(k)\Leftrightarrow X\overset{d}{=}Z_1^2+\dots+Z_k^2$, where $Z_i \overset{iid}{\sim}N(0,1)$.
"""

# ╔═╡ 46fb4c81-adf5-4b8b-83a9-01bbb37607f1
md"""
### How to generate it?
"""

# ╔═╡ ff5a1ec6-ab09-4137-84ed-068f4d9895ce
md"""
##### 자유도가 4인 카이제곱분포에서 100개의 샘플을 얻는 방법
"""

# ╔═╡ 03f2eacf-c257-48e8-867c-a1775166b98e
md"""
(방법1)
"""

# ╔═╡ e18c1f97-b009-481a-b19e-345be4146ecb
rand(Chisq(4),100)

# ╔═╡ e2db5fb2-7537-424c-8e69-e6e3247ec51a
md"""
(방법2) 정규분포 -> 카이제곱분포 
"""

# ╔═╡ 9314a572-e10a-4ed3-b882-dbbaaa17c075
[rand(Normal(0,1),4).^2 |> sum for i in 1:100]

# ╔═╡ cf98d110-505d-44af-b4c3-35e68289fe99
md"""
(방법3) 지수분포 -> 카이제곱분포
- 복습: $X,Y \overset{iid}{\sim} N(0,1) \Rightarrow R^2/2 \sim Exp(1)$, where $R^2=X^2+Y^2$.
"""

# ╔═╡ 26a5eb9e-b8ed-4580-a5e7-fd2b72f65ec3
[rand(Exponential(1))*2 + rand(Exponential(1))*2 for i in 1:100]
# rand(Exponential(1))이 R²/2 

# ╔═╡ a6c65dc1-e84c-4d35-a8f5-2d0c9db6f045
[rand(Exponential(2)) + rand(Exponential(2)) for i in 1:100]

# ╔═╡ e292a985-98e9-49e1-ac2f-4b7d696026b9
let 
	N=10000
	X1= rand(Chisq(4),N)
	X2= [rand(Normal(0,1),4).^2 |> sum for i in 1:N]
	X3= [rand(Exponential(1))*2 + rand(Exponential(1))*2 for i in 1:N]
	X4= [rand(Exponential(2)) + rand(Exponential(2)) for i in 1:N]
	histogram(X1)
	histogram!(X2)
	histogram!(X3)
	histogram!(X4)
end 

# ╔═╡ 1af05070-bb8d-45b5-a590-5900baa8b9ed
md"""
`-` 정리하면 아래와 같이 된다.
-  $Y \sim \chi^2(4)$
-  $Y \overset{d}{=} Z_1^2+Z_2^2+Z_3^2+Z_4^2$, where $Z_i \overset{iid}{\sim} N(0,1)$.
-  $Y \overset{d}{=} 2\frac{Z_1^2+Z_2^2}{2}+2\frac{Z_3^2+Z_4^2}{2}=2\frac{R_1^2}{2}+2\frac{R_2^2}{2}$, where $R_i^2/2 \overset{iid}{\sim} Exp(1)$.
-  $Y \overset{d}{=} X_1+X_2$, where $X_1,X_2 \overset{iid}{\sim} Exp(2)$.
"""

# ╔═╡ d3975d4c-9751-4a27-953c-8ee603b33eef
let
	k=6
	md"""
	`-` 자유도가 $k 인 카이제곱분포 
	- 표준정규분포 $k 개를 제곱하여 합친것과 같다. 
	- 평균이 2인 지수분포 $(k/2) 개를 합친것과 같다 (?)
	"""
end

# ╔═╡ 2030222c-b0e1-4a9d-b1fd-16ace7b47f2f
md"""
### note: 표본분산의 분포
- 이론: $X_1,\dots,X_n \overset{iid}{\sim} N(0,1) \Rightarrow \frac{(n-1)S^2}{\sigma^2}\sim \chi^2(n-1)$
- 나중에 공부합시다.. 
"""

# ╔═╡ 70532d58-2829-4e81-bdda-679c2b9dda77
md"""
### 카이제곱분포의 합
- 이론: $X \sim \chi^2(k_1),~ Y\sim \chi^2(k_2),~ X \perp Y \Rightarrow X+Y \sim \chi^2(k_1+k_2)$ 
"""

# ╔═╡ 6ba47796-f396-42b6-abcd-1df06f98da58
md"""
### 히스토그램 
"""

# ╔═╡ 83cadf65-46dd-45b7-a2e1-0b858d719676
md"자유도=$(@bind k Slider(1:200, show_value=true))"

# ╔═╡ ffe59d85-656e-4056-afa1-fb7e76b59c42
rand(Chisq(k),1000) |> histogram

# ╔═╡ 7254e500-ffaf-4280-a0cd-14318ec404aa
md"""
## 감마분포: $X \sim \Gamma(\alpha,\beta)$
"""

# ╔═╡ 80f217cc-f981-4625-a9a3-0f8450ab0d13
md"""
### 감마분포 요약
- X의의미: 서로 독립인 지수분포를 $\alpha$개 합친 것, 시간1에 평균적으로 $\lambda$번 발생하는 사건이 있을때 $\alpha$번째 사건이 발생할때까지 걸리는 시간.
- X의범위: $x \in (0,\infty)$
- 파라메터의 의미: $\alpha=$ 지수분포를 더한 횟수(의 확장버전), $\beta=\frac{1}{\lambda}=$ 지수분포의 평균
- 파라메터의 범위: $\alpha>0$, $\beta>0$.
- pdf: $\frac{1}{\Gamma(\alpha)\beta^{\alpha}}x^{\alpha-1}e^{-x/\beta}$
- mgf: 
- E(X): $\alpha\beta$
- V(X): $\alpha\beta^2$
"""

# ╔═╡ 04d93dc6-43d3-4081-bca6-9c08f5ee2653
md"""
### 대의적 정의 ($\alpha$가 자연수일경우)
`-` $X \sim \Gamma(\alpha,\beta) \Leftrightarrow X \overset{d}{=} Z_1+\dots+Z_\alpha$, where $Z_i \overset{iid}{\sim} Exp(\frac{1}{\beta})$
"""

# ╔═╡ 36e5c715-ec53-4586-b006-7e51ae83255e
md"""
-  $Exp(\frac{1}{\beta})$는 평균이 $\beta$인 지수분포 
"""

# ╔═╡ a3894a73-7fee-43ec-a06a-ff0142ccacc6
md"""
### how to generate it?
"""

# ╔═╡ 7130347c-0a73-4f7d-8c64-2164eb65b37d
md"""
##### $\Gamma(3,2)$를 1000개 생성하라. 
"""

# ╔═╡ e4e7bbcf-7e48-4020-9e19-4cbd2d7d36e3
md""" 
α = $(@bind α Slider(1:100, show_value=true, default=3))

β = $(@bind β Slider(0.1:0.1:100, show_value=true, default=2))
"""

# ╔═╡ 5662230c-64a0-4835-b8a6-7dcf3fe6ae55
md"""
(방법1)
"""

# ╔═╡ 42f1313b-528d-4aec-85bb-86ff6bbb288a
rand(Gamma(α,β),1000)

# ╔═╡ 4a21db7e-2e84-438f-a8ba-a20360628777
md"""
(방법2) 지수분포 -> 감마분포 
"""

# ╔═╡ 7b37ed60-e13a-41a7-9efb-1c26c9017977
[rand(Exponential(β),α) |> sum for i in 1:1000]

# ╔═╡ 67508ce7-77d8-4ce0-9a0b-0312ffca76ec
md"""
(방법3) 표준정규분포 이용
"""

# ╔═╡ 2bdeef3c-ad0c-4dcd-8d45-a8d0e88a743f
[rand(Normal(0,1),2α).^2 |> sum for i in 1:1000]
# 방법3이 가능한 이유는 2α가 현재 자연수이고 β=2이기 때문 

# ╔═╡ f6d6ca0c-ad73-400a-93d7-086453df8470
md"""
(방법4) 카이제곱분포를 이용
"""

# ╔═╡ bb13a418-685c-41c1-972d-cf2b2a978ec8
rand(Chisq(2α),1000)
# 방법4이 가능한 이유는 2α가 현재 자연수이고 β=2이기 때문 

# ╔═╡ 2c195c00-bee6-49d5-9f51-459f37bcebcb
md"""
`-` 확인
"""

# ╔═╡ 195ef5c7-a4fa-4176-a0a7-846ee354ce7f
let
	N=10000
	X1= rand(Gamma(α,β),N)
	X2= [rand(Exponential(β),α) |> sum for i in 1:N]
	X3= [rand(Normal(0,1),2α).^2 |> sum for i in 1:N]
	X4= rand(Chisq(2α),N)
	histogram([X1,X2,X3,X4])
end 

# ╔═╡ c01d6d44-9310-46c1-b040-b22bc19b936b
md"""
### 감마분포와 카이제곱분포의 관계
"""

# ╔═╡ cc16de8e-d845-448d-8dbe-c1e8988d48df
md"""
`-` 이론: $X \overset{d}{=} Y$, where $X\sim \chi^2(k)$ and $Y\sim \Gamma(\frac{k}{2},2)$.
"""

# ╔═╡ 906136a9-095c-4ae5-a75d-290f14e2f1a6
let
	k=3 
	histogram(rand(Chisq(k),10000))
	histogram!(rand(Gamma(k/2,2),10000))
end 

# ╔═╡ dea33225-75fd-448f-9854-aef9d2973a47
md"""
### 히스토그램 
"""

# ╔═╡ 52606dc6-6e8b-4ab3-88e5-ef0f3bb73fbd
md"""
α = $(@bind shape Slider(0.1:0.1:30))
β = $(@bind scale Slider(0.1:0.1:100))
"""

# ╔═╡ 4415ff22-04a9-48d5-95c8-333921469951
histogram(rand(Gamma(shape,scale),10000000), title=md"Gamma($shape, $scale)")

# ╔═╡ f8f2dd38-772f-4afb-9fd6-847f43bce147
md"""
### 척도모수 
"""

# ╔═╡ c2e4e893-4590-49ad-be6b-39731fa1bdd9
md"""
`-` 감마분포도 척도모수를 가짐 
"""

# ╔═╡ 0456d788-ffdb-4457-be96-52a2f4704084
begin
	histogram(rand(Gamma(50,6),10000))
	histogram!(rand(Gamma(50,2)*3,10000))
end 

# ╔═╡ 0cf5faae-20c8-441c-a98f-3af7789e5367
md"""
### 감마분포의 합 
"""

# ╔═╡ a37b1e80-943b-4db1-9231-0fb0dcd161ce
md"""
`-` 이론: $X \sim \Gamma(\alpha_1,\beta), ~ Y\sim \Gamma(\alpha_2,\beta),~ X \perp Y \Rightarrow X+Y \sim \Gamma(\alpha_1+\alpha_2,\beta)$
"""

# ╔═╡ a4f3b441-91b7-445f-b516-2a20424092a1
begin
	histogram(rand(Gamma(25,6),10000))
	histogram!(rand(Gamma(6,6),10000)+rand(Gamma(25-6,6),10000))
end

# ╔═╡ c296c866-08b9-4a59-bf6b-e5f8c162c2b9
md"""
## 지수분포의 다양한 표현
"""

# ╔═╡ c1fe1819-20a1-47a8-b96e-bbbd6fec9b00
md"""
`-` 지수분포는 따라서 다양하게 표현된다.
- 경우1: $X\sim Exp(\theta)$, $f(x)=\frac{1}{\theta}e^{-\frac{x}{\theta}}$.
- 경우2: $X\sim Exp(1/\lambda)$, $f(x)=\lambda e^{-\lambda x}$.
- 경우3: $X\sim Exp(\lambda)$, $f(x)=\lambda e^{-\lambda x}$.

`-` 기억할것 
- 경우1: 지수분포의 모수는 평균, 지수분포의 파라메터는 $\theta=\frac{1}{\lambda}$으로 정의하여 새롭게 사용 
- 경우2: 지수분포의 모수는 평균, 지수분포의 파라메터는 포아송의 $\lambda$를 재활용. 
- 경우3: 지수분포의 모수는 평균의 역수, 지수분포의 파라메터는 포아송의 $\lambda$를 재활용. 

`-` 노테이션의 숨은의도들 (제 생각)
- 경우2: 포아송분포의 파라메터도 그대로 쓰고싶고, "지수분포의 모수 = 지수분포의 평균" 과 같이 만들고 싶음 
- 경우1: 경우2에서는 $X \sim Exp(1/\lambda)$로 표현되어서 모수가 역수로 되어있어 헷갈림. 그냥 포아송분포의 $\lambda$를 버리는 편이 좋겠음. 지수분포의 평균을 의미하는 $\theta$를 새롭게 정의하고 이 $\theta$를 중심으로 pdf를 만듬 
- 경우3: 경우2에서는 $X \sim Exp(1/\lambda)$로 표현되어서 모수가 역수로 되어있어 헷갈림. 그냥 모수는 그대로 $Exp(\lambda)$를 쓰고 지수분포의 평균을 모수의 역수로 정의하는 것이 좋겠음. 

`-` 아무튼 여러가지 방식으로 표현합니다..
- 경우1: 줄리아, 파이썬, 위키
- 경우2: 제가 공부한 교재? 요즘은 이렇게 잘 안쓰는 것 같음 
- 경우3: R, 위키

`-` 평균이 $314$인 지수분포의 pdf는 $f(x)=\frac{1}{314}e^{-x/314}$ 이다. <- 이렇게 외웠음
"""

# ╔═╡ ca700e4d-5810-4c66-8347-2fe2463ccc59
md"""
## 숙제
`-` 평균이 2인 지수분포를 이용하여 자유도가 20인 카이제곱분포를 생성하여라.
"""

# ╔═╡ f2563532-8315-4439-9aed-6a62753eeb73
let 
	N=10000
	[sum(rand(Exponential(2),10)) for i in 1:N] |> histogram
	rand(Chisq(20),N) |> histogram!
	#[sum(rand(Chisq(2),10)) for i in 1:N] |> histogram!
end 

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Distributions = "~0.25.53"
Plots = "~1.27.5"
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
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

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
git-tree-sha1 = "621f4f3b4977325b9128d5fae7a8b4829a0c2222"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.4"

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
git-tree-sha1 = "88ee01b02fba3c771ac4dce0dfc4ecf0cb6fb772"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.27.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "670e559e5c8e191ded66fa9ea89c97f10376bb4c"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.38"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "d3538e7f8a790dc8903519090857ef8e1283eecd"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.5"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

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
git-tree-sha1 = "4f6ec5d99a28e1a749559ef7dd518663c5eca3d5"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "8d7530a38dbd2c397be7ddd01a424e4f411dcc41"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.2"

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
# ╟─b8113ace-bf87-11ec-332c-d57eb9d49518
# ╟─1cfafdd4-86ff-4e2f-aa50-e9e78f4f4090
# ╟─e72224fa-ebba-48e4-b702-b2d0329610a6
# ╠═9d2a3140-7c5e-4439-93e0-1e02e8417daf
# ╠═7fd8042e-9af6-4f6c-b25c-ab8ed5f0af93
# ╠═aa3f4785-8e3f-4841-b3cd-fb1077792992
# ╟─fe50793a-c255-4cc6-81fe-21d729262786
# ╟─823932ac-0546-427f-8d1a-bb1a452ef29c
# ╟─019f20b3-7920-40e9-b470-4f6392b1c4ac
# ╟─930005f3-d5dc-4642-ab96-04f9c74a1c8f
# ╠═0fd62af9-1519-4a64-a860-99e862087320
# ╠═410876fb-e8d3-413b-b278-9c0d7e76b9a4
# ╟─01e17c60-fd70-4767-b51b-cfd7d542df9f
# ╟─90306510-4d75-4189-89f1-97d3a1f83f21
# ╠═239b3846-48df-40b1-b223-fbffa87cd5a4
# ╟─076614a8-3c5a-47cc-bcc8-c4074e9776fe
# ╠═528f54ca-8c1a-4e99-84bb-e1c91004057f
# ╠═f7274989-9cd5-48dd-b501-91719cc07703
# ╠═170bf0c7-3055-4c7c-95e2-e2d3d9c3244c
# ╟─69d57189-a1d8-4574-9714-d5524ff94c55
# ╠═1912ea52-b536-4200-91a1-d9580e635b93
# ╟─c8a2897c-1381-4828-8d87-d9e3140788e6
# ╠═1d5ced10-a8b3-4e2c-9f0f-e2b7d274c8c3
# ╟─3dd308dd-3f64-4a1c-a0e5-96d88d1a68c7
# ╟─abf6b2ce-b5b0-4921-8edd-bda6b38ff364
# ╟─67d81f4b-4b16-4b0b-aacb-1d38a6e1dbc5
# ╟─9061c67f-5039-4ee2-a749-d51aa8034320
# ╟─7fa5ef2f-6aa2-4c00-bf4b-fc47c409dd10
# ╟─a46d857a-6328-43c2-8aff-7155462c2119
# ╟─27a1f077-e674-40ab-889c-06663a133c8a
# ╟─46fb4c81-adf5-4b8b-83a9-01bbb37607f1
# ╟─ff5a1ec6-ab09-4137-84ed-068f4d9895ce
# ╟─03f2eacf-c257-48e8-867c-a1775166b98e
# ╠═e18c1f97-b009-481a-b19e-345be4146ecb
# ╟─e2db5fb2-7537-424c-8e69-e6e3247ec51a
# ╠═9314a572-e10a-4ed3-b882-dbbaaa17c075
# ╟─cf98d110-505d-44af-b4c3-35e68289fe99
# ╠═26a5eb9e-b8ed-4580-a5e7-fd2b72f65ec3
# ╠═a6c65dc1-e84c-4d35-a8f5-2d0c9db6f045
# ╠═e292a985-98e9-49e1-ac2f-4b7d696026b9
# ╟─1af05070-bb8d-45b5-a590-5900baa8b9ed
# ╟─d3975d4c-9751-4a27-953c-8ee603b33eef
# ╟─2030222c-b0e1-4a9d-b1fd-16ace7b47f2f
# ╟─70532d58-2829-4e81-bdda-679c2b9dda77
# ╟─6ba47796-f396-42b6-abcd-1df06f98da58
# ╠═83cadf65-46dd-45b7-a2e1-0b858d719676
# ╠═ffe59d85-656e-4056-afa1-fb7e76b59c42
# ╟─7254e500-ffaf-4280-a0cd-14318ec404aa
# ╟─80f217cc-f981-4625-a9a3-0f8450ab0d13
# ╟─04d93dc6-43d3-4081-bca6-9c08f5ee2653
# ╟─36e5c715-ec53-4586-b006-7e51ae83255e
# ╟─a3894a73-7fee-43ec-a06a-ff0142ccacc6
# ╟─7130347c-0a73-4f7d-8c64-2164eb65b37d
# ╟─e4e7bbcf-7e48-4020-9e19-4cbd2d7d36e3
# ╟─5662230c-64a0-4835-b8a6-7dcf3fe6ae55
# ╠═42f1313b-528d-4aec-85bb-86ff6bbb288a
# ╟─4a21db7e-2e84-438f-a8ba-a20360628777
# ╠═7b37ed60-e13a-41a7-9efb-1c26c9017977
# ╟─67508ce7-77d8-4ce0-9a0b-0312ffca76ec
# ╠═2bdeef3c-ad0c-4dcd-8d45-a8d0e88a743f
# ╟─f6d6ca0c-ad73-400a-93d7-086453df8470
# ╠═bb13a418-685c-41c1-972d-cf2b2a978ec8
# ╟─2c195c00-bee6-49d5-9f51-459f37bcebcb
# ╠═195ef5c7-a4fa-4176-a0a7-846ee354ce7f
# ╟─c01d6d44-9310-46c1-b040-b22bc19b936b
# ╟─cc16de8e-d845-448d-8dbe-c1e8988d48df
# ╠═906136a9-095c-4ae5-a75d-290f14e2f1a6
# ╟─dea33225-75fd-448f-9854-aef9d2973a47
# ╠═52606dc6-6e8b-4ab3-88e5-ef0f3bb73fbd
# ╠═4415ff22-04a9-48d5-95c8-333921469951
# ╠═f8f2dd38-772f-4afb-9fd6-847f43bce147
# ╠═c2e4e893-4590-49ad-be6b-39731fa1bdd9
# ╠═0456d788-ffdb-4457-be96-52a2f4704084
# ╟─0cf5faae-20c8-441c-a98f-3af7789e5367
# ╟─a37b1e80-943b-4db1-9231-0fb0dcd161ce
# ╠═a4f3b441-91b7-445f-b516-2a20424092a1
# ╟─c296c866-08b9-4a59-bf6b-e5f8c162c2b9
# ╟─c1fe1819-20a1-47a8-b96e-bbbd6fec9b00
# ╟─ca700e4d-5810-4c66-8347-2fe2463ccc59
# ╠═f2563532-8315-4439-9aed-6a62753eeb73
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
