### A Pluto.jl notebook ###
# v0.19.40

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

# ╔═╡ 521890de-ab23-11ec-0c2f-2dcaee6dc1bc
using Plots, Distributions, PlutoUI,Random

# ╔═╡ dbf8cdd6-7815-49f8-9dd6-5987000792ce
md"""
# 03wk-1: 지수분포, 박스뮬러변환
"""

# ╔═╡ b7a1c25e-1aeb-4a57-80b5-a9ea4c5c0530
md"""
## 1. 강의영상
"""

# ╔═╡ 20697867-c74a-4485-bb97-c6a31060669a
# html"""
# <div style="display: flex; justify-content: center;">
# <div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
# <iframe src=
# "
# https://www.youtube.com/embed/playlist?list=PLQqh36zP38-xIz4hvPmHZz0g8g04MJ0xZ
# "
# width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
# """

# ╔═╡ 64e29d20-aaf5-4fa0-ad03-b283dac52dce
md"""
## 2. Imports
"""

# ╔═╡ 4e838664-c6ea-4421-99e6-c585fd21cb12
Plots.plotly()

# ╔═╡ 61967452-a8e0-46a8-85cb-c0944451c5d2
PlutoUI.TableOfContents(aside=true)

# ╔═╡ 29a2e115-37a7-4f67-a040-a4c1648fdfa6
md"""
## 3. 지수분포 ($X \sim Exp(\lambda)$)
"""

# ╔═╡ a185e773-f51e-449b-8d40-de819bd4badf
md"""
### A. 기본내용
"""

# ╔═╡ 78ebbf73-8663-490b-8c0f-58bba2fe2ead
md"""
-- 간단한 요약 
- X의의미: 시간1에 평균적으로 λ번 발생하는 시건이 있을때 첫 번째 이벤트가 발생할때 까지 걸리는 시간. 
- X의범위: 시간은 양수이므로 X ≥ 0
- 파라메터의 의미: (1) λ = 시간1에 평균적으로 발생하는 이벤트의 수 (2) 1/λ = θ = 한번의 이벤트가 발생할때까지 평균적으로 걸리는 시간 
- 파라메터의 범위: λ>0, θ>0. 
- pdf: $f(x)=\lambda e^{-\lambda x}=\frac{1}{\theta}e^{-\frac{1}{\theta}x}$ 
- mgf: 
- cdf: $F(x)=1-e^{-\lambda x}=1-e^{-\frac{1}{\theta}x}$
- E(X): $\frac{1}{\lambda}=\theta$
- V(X): $\frac{1}{\lambda^2}=\theta^2$
"""

# ╔═╡ 36af80d2-bff4-426e-8b44-8f407638941f
md"""
-- 지수분포의 다양한 표현

- 경우1 (위키, 줄리아, 파이썬): $X \sim Exp(\theta), f(x)=\frac{1}{\theta}e^{-\frac{x}{\theta}}$
- 경우2 (위키, R): $X \sim Exp(\lambda), f(x)=\lambda e^{-\lambda x}$
- 경우3: $X \sim Exp(1/\lambda), f(x)=\lambda e^{-\lambda x}$
"""

# ╔═╡ 770f89b0-25ac-42ed-a320-74809e368b44
md"""
!!! info "약속"
	앞으로는 경우1,3을 혼용하여 사용한다. 경우2는 사용하지 않는다 (줄리아 문법상 충돌). 헷갈림을 방지하기 위해서 지수분포를 의미할땐 평균이 $\theta$ 인 지수분포라고 표현하자. 아래와 같이 pdf와 cdf도 기억하자. 
	- pdf: ``f(x)=\frac{1}{\theta}e^{-\frac{1}{\theta}x}=\lambda e^{-\lambda x}``
	- cdf: ``F(x)=1-e^{-\frac{1}{\theta}x}=1-e^{-\lambda x}``
	줄리아에서의 코드는 아래와 같이 기억하자. 
	- 줄리아에서: `Exponetial(θ)`, `Exponetial(1/λ)`
"""

# ╔═╡ 68f0c9b7-25e9-4807-a543-0c7a9a554898
md"""
### B. 모수 $\to$ 히스토그램
"""

# ╔═╡ 51231561-64d1-4545-93d7-980bf77847ae
# 슬라이더 바 생성

# ╔═╡ dea01e57-ca18-4a54-9631-8ca3046f4bd8
md"""
*Fig -- 평균이 ``\theta`` 인 지수분포 (왼쪽: 히스토그램, 오른쪽: pdf)*
"""

# ╔═╡ 39235f8a-4e2d-4b26-82f4-d95a5c2167ed
# 

# ╔═╡ 885941da-8604-4af0-8b69-6d564348f116
md"""
### C. 난수생성 테크닉
"""

# ╔═╡ cb63ae3e-7de4-4eac-aca3-55561864956d
md"""
*평균이 1/4인 지수분포에서 100개의 샘플을 뽑는 방법*
"""

# ╔═╡ 2134273b-fca8-4238-b906-f05948fae7d2
md"""
(방법1)
"""

# ╔═╡ 1da1e41d-73a6-43f4-9d69-3e28b83da8e0
rand(Exponential(1/4),100)

# ╔═╡ caa4ea34-e2ae-455a-88a7-6324ff2cc8ef
md"""
(방법2) 포아송 -> 지수분포 (X), 포아송프로세스 -> 지수분포 (O)

-- (틀린)아이디어: ``X``를 포아송분포에서 추출한뒤에 ``\frac{1}{X}``을 하면 지수분포 나오는거 아니야?

- 맥도날드에 시간1당 4명씩 평균적으로 방문한다. 1명 방문하는데에는 평균적으로 시간이 0.25 걸린다고 볼 수 있음. 
- 따라서 언뜻생각하면 포아송과 지수분포는 역의 관계라서 포아송분포를 만들고 역수를 취하면 지수분포를 쉽게 만들 수 있을 것 같다.
"""

# ╔═╡ 7620f927-e031-44f5-b22e-bfaa26bdad84
md"""
-- (틀린)아이디어를 시도하기 위해서 일단 포아송을 뽑았는데..
"""

# ╔═╡ eb7572bd-63eb-4792-b982-9f320358a99c
rand(Poisson(4),100)

# ╔═╡ 1c902aac-d633-4cde-95e1-62e42e0ffe87
md"""
- 0이 나온다?
- 생각해보니까 0이 없다고 쳐도 나올 수 있는 값은 1, 1/2, 1/3, 1/4, ... 따위임 (애초에 틀린 접근)
"""

# ╔═╡ 1e417760-f0c7-4dad-8ffd-003a1ea24f06
md"""
-- 새로운 아이디어: 포아송프로세스 셋팅 $\to$ 지수분포

- 극한의 베르누이로 포아송을 만들때, **"몇번 성공했는지"** 관심을 가지고 카운팅 했다. 
- 이를 조금 응용해서 **"첫 성공까지 몇번의 시도를 해야하는지"** 카운팅을 한다고 생각하면 시간계산이 가능할것 같다. 
- 결국 "포아송분포 -> 지수분포"로 추출하는 것이 아니라 "포아송프로세스 -> 지수분포"와 같은 방식으로 추출해야 한다. 

"""

# ╔═╡ 69324599-c085-45f9-8e6b-359c2d1a163c
md"""
(예비학습) 기하분포 
"""

# ╔═╡ 5d39448b-2e68-475d-a6d4-a06c3b6b347a
# 성공할때까지 시도하는 함수: 성공확률 -> 1회성공까지 시도한 횟수

# ╔═╡ f3ee3d9c-4163-469a-b96f-8536a413025c
md"""
*Fig -- (a) 줄리아에서 기하분포 추출 (b) `mygeo`에서 기하분포 추출*
"""

# ╔═╡ fae9c055-f4dd-4d6e-9dc7-ddd6a74fd905


# ╔═╡ aef02ab9-4e12-4030-b428-f4aa36f65ba6
md"""
(풀이시작)
"""

# ╔═╡ 6924193e-ea15-42e9-a702-8c6b0a224f35
md"""
*Fig -- 포아송프로세스를 응용해서 지수분포를 생성*
"""

# ╔═╡ a27efaae-2793-4d8a-8dd6-4900caa80a4a


# ╔═╡ fc1acb74-f856-47c7-b9e1-8d635f0e35ca
md"""
(방법3) inverse cdf method 
- 이론적인 pdf를 알고 있다는 전제가 필요함. 
- 자세하게 살펴보자. 
"""

# ╔═╡ c2c198a7-4597-49e2-98d7-8776603db7fa
md"""
*Inverse cdf method를 활용하여 지수분포에서 샘플추출*
"""

# ╔═╡ cf6a955e-d684-4a16-bc0f-14409f83536a
md"""
`-` 아래와 같은 2개의 지수분포의 pdf를 고려하자.

$$f(x)=e^{-x}$$ 

$$g(x)=\frac{1}{5}e^{-\frac{1}{5}x}$$
"""

# ╔═╡ ed777a37-6351-41a1-aed3-86412c9edaac
md"""
`-` 각각의 pdf를 그려보면 아래와 같다. 
"""

# ╔═╡ 8ab6d3db-1037-4afb-9d30-d16675a1bffd
md"""
*Fig -- 평균이 각각 1,5인 지수분포의 pdf*
"""

# ╔═╡ db26ed1c-b01a-4f06-886d-4079cf2c139f


# ╔═╡ 7889eeaa-a7b5-4cdd-8b26-c0de2d2981e0
md"""
`-` 이번에는 각각의 cdf를 그려보자. 

-  $F(x)=\int_0^x f(\tau)d\tau=\int_0^x e^{-\tau} d\tau = \left[-e^{-\tau}\right]_0^x=1-e^{-x}$
-  $G(x)=\int_0^x g(\tau)d\tau=\int_0^x \frac{1}{5}e^{-\tau/5} d\tau = \left[-e^{-\tau/5}\right]_0^x=1-e^{-x/5}$

"""

# ╔═╡ 0e381090-dac3-48fe-adbb-f25e630d25bf
md"""
*Fig -- 평균이 각각 1, 5인 지수분포의 cdf*
"""

# ╔═╡ 23abf65d-276d-4a5d-ae70-b41255b35b48


# ╔═╡ ac1c2222-04b5-4af0-a55c-723b1ad57dec
md"""
`-` cdf 해석 
- 왼쪽(평균이1인지수분포): 5정도면 거의 cdf의 값이 1에 가까워짐 
- 오른쪽(평균이5인지수분포): 5정도에서 값이 거의 0.63정도임 $\to$ 100번뽑으면 5보다 작은게 63개정도.. 

`-` cdf의 y축에서 랜덤변수를 발생시킨다음 $\rightarrow \downarrow$ 와 같이 이동하여 $x$축에 내린다고 생각해보자. 
- 왼쪽: 대부분 5이하에 떨어짐 
- 오른쪽: 약 63% 정도만 5이하에 떨어짐.
"""

# ╔═╡ e4a7a78e-e2c5-4332-ba6e-efe0a87dd4c8
md"""
*Fig -- Inverse cdf를 이해하기 위한 아이디어*
"""

# ╔═╡ 9c2f4080-9dd5-4f85-83f7-42bf3c719e6e


# ╔═╡ e234260f-b0bd-4f0e-a5cd-bdc9afad041f
md"""
- 빨간색: $X \sim$ 균등분포
- 초록색: $X \sim$ 지수분포???
"""

# ╔═╡ cd6507ca-e30c-4c2c-8d4a-11b714db833e
md"""
*Fig -- Inverse cdf를 이용하여 지수분포 샘플링*
"""

# ╔═╡ b5093fca-2946-41b2-a023-d4abc17c453c


# ╔═╡ d7823cc8-231f-41b4-8a33-175cd42cb2b5
md"""
### D. 무기억성 
"""

# ╔═╡ 35760f05-e335-46cb-baea-83763f3f20d5
md"""
!!! info "이론: 지수분포의 무기억성"
	평균적으로 10분마다 한번씩 이벤트가 발생하는 어떠한 사건이 있다고 하자. (평균이 10분이 지수분포를 따르겠지?) 이때 9분59초까지 이벤트가 발생하지 않았다고 가정하자. 그렇다고 하여도 "곧 약속의 시간 10분이 되어서 이벤트가 발생한다" 는 보장은 없다. 왜냐하면 

	-  $X \sim Exp(\lambda) \Rightarrow$ 모든 $t,s>0$에 대하여 $P(X>t)=P(X>t+s|X>s)$

	이기 때문이다. 
"""

# ╔═╡ fbb0395a-9ea9-4eb6-a439-51aa13024a31
md"""
-- 개념: 

- 이해를 위해서 $t=1,s=9$ 대입 => $P(X>1)=P(X>10 | X>9)$
- 좌변: 시간을 1 기다려서 이벤트가 발생안할 확률 
- 우변: 시간을 9 기다렸는데 이벤트가 발생안했음 -> 시간 10을 기다려서 이벤트가 발생안할 확률 
- 예를들어서 $\lambda=0.1$ 이라면 한번 이벤트 발생하는데 평균 시간10이 걸린다는 의미임. => (1) 좌변은 이제 시간1기다림 (2) 우변은 시간 9를 기다림. 곧 "약속된" 시간 10이 완성됨 => 우변이 더 확률이 크지 않을까? => 아니라는거에요!
"""

# ╔═╡ 68dcc6c9-7c52-4d16-b5b0-f0aef055c01a
md"""
-- 이해: 지수분포의 근본? 포아송 프로세스 
- 엄청 짧은 시간 
- 엄청 작은 확률 
- 엄청 많은 베르누이 시행이 "독립적"으로 수행됨 -> 지금까지 실패했다고 해서 이후에 성공확률이 높아지는건 아님. 
- 우변: 이미 시간9동안 무수히 많은 독립적인 베르누이 시행을 놓친상태임. 그 이후의 시행은 모두 독립이므로 좌변의 확률보다 더 크다고 볼 수 없음. 
"""

# ╔═╡ 45a26fc2-3c96-4130-83b2-e75c39846264
md"""
-- 간단한 실습1
"""

# ╔═╡ 86f22bd5-fb69-46cc-b53e-e1e04e99b5e2
md"t= $@bind t Slider(0.01:0.01:5,show_value=true,default=1.0)" 

# ╔═╡ f455bad9-640f-4f6d-ab08-04d432681de9
md"s= $@bind s Slider(0.01:0.01:5,show_value=true,default=2.0)" # 이미 기다린시간

# ╔═╡ 6de8f87f-7f99-480a-b5f6-30f0a152f5a1
let
	N = 5000
	X = rand(Exponential(1),N) 
	println("P(X>t) = $(sum(X .> t) / N)")
	println("P(X>t+s|X>s) = $(sum(X .> t+s) / sum(X .> s))")
end

# ╔═╡ 67edeecf-a93a-488d-844b-450866d9b43f
md"""
-- 간단한 실습2 

- 무기억성 = 과거는 중요하지 않음! 
-  $P(X>1)=P(X>2|X>1)=P(X>3|X>2)=...$ 
"""

# ╔═╡ 5452a29a-2e5a-47c2-804f-7a5cad2d943f
let 
	N = 5000
	X = rand(Exponential(1),N)
	println("P(X>1) = $(sum(X.>1)/N)")
	println("P(X>2|X>1) = $(sum(X.>2)/sum(X.>1))")
	println("P(X>3|X>2) = $(sum(X.>3)/sum(X.>2))")
end

# ╔═╡ 72edfa85-7ce4-49ea-a11a-49d714536106
md"""
### E. 척도모수 ($\star$)
"""

# ╔═╡ b0fefc2a-a6b6-4615-8a94-98a328a09042
md"""
!!! info "이론: 지수분포는 척도모수를 가진다."
	평균이 1인 지수분포를 $a$배하면 평균이 $a$인 지수분포를 따른다. 즉 아래가 성립한다. 
	-  $X \sim Exp(θ) \Rightarrow aX \sim Exp(a\theta)$
    이를 이용하면 평균이 2인 지수분포는 평균이 1인 지수분포를 만든뒤에 2를 곱하면 된다.
"""

# ╔═╡ b66b5d63-9e49-4456-88c3-d47c48870fe9
md"""
*Fig -- 지수분포가 척도모수를 가진다는 사실을 알면, 너무 편하게 샘플들을 추출할 수 있음.* 
"""

# ╔═╡ ecbd2a2b-e527-436a-95a1-e62d85cd91b0
let 
	N = 1000
	rand(Exponential(12),N) |> histogram
	rand(Exponential(1),N).*12 |> histogram!
	rand(Exponential(2),N).*6 |> histogram!
	rand(Exponential(3),N).*4 |> histogram!
	rand(Exponential(4),N).*3 |> histogram!
	rand(Exponential(6),N).*2 |> histogram!
end 

# ╔═╡ 713adf30-a586-47c7-832c-5b6588203862
md"""
## 4. 박스뮬러변환
"""

# ╔═╡ e88df4b7-5352-49e4-be65-6be16f23f003
md"""
### A. 정규분포와 지수분포의 관계
"""

# ╔═╡ 5d024a81-0244-4348-b6c6-e5e9f70545bd
md"""
!!! info "이변량 정규분포와 지수분포의 관계"
	``\begin{bmatrix} X \\ Y \end{bmatrix} \sim N({\bf 0},{\bf I})`` 일때 반지름제곱 ``R^2=X^2+Y^2``은 평균이 2인 지수분포를 따른다. 즉 
	
	- ``\begin{bmatrix} X \\ Y \end{bmatrix} \sim N({\bf 0},{\bf I}) \quad \Rightarrow \quad \begin{cases} X^2+Y^2 \sim Exp(2) \\ \Theta \sim U(0,2\pi) \end{cases} ``
	
	이는 서로 독립인 2개의 표준정규분포로 지수분포를 만들 수 있다는 사실을 의미한다. 또한 역으로 아래도 성립한다. 
	
	- ``\begin{cases} R^2/2 \sim Exp(1) \\ \Theta \sim U(0,2\pi) \end{cases} \Rightarrow \begin{bmatrix} R\cos \Theta \\ R \sin \Theta \end{bmatrix} \sim N\left (\begin{bmatrix} 0 \\ 0 \end{bmatrix}, \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}\right)``
	
	이는 지수분포와 균등분포로 정규분포를 만들 수 있다는 사실을 의미한다.
"""

# ╔═╡ 6cf0f4d7-d8f4-4478-8b30-55fbeaeee150
md"""
-- 의미?
"""

# ╔═╡ bb42b213-ac3d-4156-9d12-445e0702968b
@bind i Slider(1:5000,show_value=true)

# ╔═╡ 413f316d-13a4-4bda-9867-b2439be20f5c
md"""
*Fig -- 정규분포와 지수분포의 관계 (이 그림 외우세요!)*
"""

# ╔═╡ 1f99cbe2-e472-4916-80bc-ca290ca17f8b
let 
	Random.seed!(43052)
	X = rand(Normal(0,1),5000)
	Y = rand(Normal(0,1),5000)
	xi,yi = X[i],Y[i]
	R² = xi^2 + yi^2
	p1 = scatter(X,Y,alpha=0.1,label= "(X,Y), X~정규, Y~정규")
	plot!([0,xi],[0,yi],linewidth=5,label="R² = X² + Y²")
	title!("R² = $(R²)")
	p2 = (X.^2 + Y.^2) |> x-> histogram(x,alpha=0.5, label="왼쪽그림의 반지름")
	rand(Exponential(2),5000) |> x-> histogram!(x,alpha=0.5,label="평균2인 지수분포")
	plot(p1,p2)
end

# ╔═╡ 1cbdb4a1-0ddb-41e3-8130-ee79005cc193
md"""
### B. 박스뮬러변환
"""

# ╔═╡ 055bcec8-5d53-4888-bf46-50eef3a3e821
md"""
*Fig: 지수분포 $\to$ 정규분포*
"""

# ╔═╡ c565e07f-7bfb-4dfa-9b90-90171d6b5591
let 
	N = 10000
	R = .√(2*rand(Exponential(1),N))
	Θ = rand(N).*2π
	T = (R,Θ) -> (R*cos(Θ), R*sin(Θ))
	XY = T.(R,Θ) 
	X = [XY[i][1] for i in 1:N]
	Y = [XY[i][2] for i in 1:N]
	#scatter(X,Y)
	p1=histogram(X)
	p2=randn(N) |> histogram
	plot(p1,p2,layout=(2,1),xlim=(-5,5))
end 

# ╔═╡ d31dbd9c-77c6-4c18-802a-c0967ed71cab
md"""
-- inverse cdf 기법과 결합하면?
"""

# ╔═╡ f23c6456-e364-4658-a71f-6be8388209d0
md"""
`-` inverse cdf 기법과 합치면 아래와 같이 정리가능하다. 

$\begin{cases}
X=\sqrt{-2\log(1-U_1)} \cos(2\pi U_2) \\ 
Y=\sqrt{-2\log(1-U_1)} \sin(2\pi U_2) 
\end{cases},~ U_1,U_2 \overset{iid}{\sim} U(0,1)$
"""

# ╔═╡ 859d98fe-e49e-45d3-9a6f-a6354413ab91
md"""
*Fig: 균등분포 $\to$ 정규분포 (박스뮬러변환)*
"""

# ╔═╡ b944f877-0984-425b-b35c-09f1d71ffb05
let 
	N = 10000
	U1 = rand(N)
	U2 = rand(N)
	X = @. √(-2log(1-U1))*cos(2π*U2)
	p1 = histogram(X) 
	p2 = histogram(randn(N))
	plot(p1,p2,layout=(2,1),xlim=(-5,5))
end 

# ╔═╡ 6e41df13-d03c-4def-a5b1-2166c945a626
md"""
## A1. Inverse CDF의 이론적 근거
"""

# ╔═╡ f290b98a-31bb-4526-9b56-4543d2108400
md"""
!!! info "이론: Inverse CDF의 이론적 근거"
	어떠한 확률변수 $X$를 뽑고 싶다면 균등분포에서 하나의 난수를 생성하고 cdf의 역함수를 취하면 된다. 즉 아래가 성립한다. 
	-  $F^{-1}(U)\overset{d}{=} X$ 
	이때 $X \sim F$ 이고 $U\sim U(0,1)$ 이다. 
"""

# ╔═╡ 6bf2ec59-4457-4daf-b143-24b0f8ee4e0e
md"""
(증명) 이 수업에서는 증명의 편의성을 위하여 $X$가 연속형이고 $F$가 순증가함수임을 가정한다. 하지만 일반적인 경우에도 성립한다.  
"""

# ╔═╡ 0e963046-9cea-43cc-a4f4-491b59d27845
md"""
**1단계:** 먼저 $F(X)\overset{d}{=}U$임을 보이자. 
"""

# ╔═╡ b9a9e616-51d2-409e-863e-2aa9b132a75f
md"""
*Fig: 임의의 확률변수 ``X``에 대한 ``F(X)``의 히스토그램, 단 ``F(x)``는 ``X``의 cdf*
"""

# ╔═╡ de13ff8f-98ec-450c-a609-79b9683afd03
let 
	rand(Exponential(1),10000) .|> (x -> 1-exp(-x)) |> histogram
end

# ╔═╡ cd1a64c8-51f7-423e-8fa2-d67351771482
md"""
-- 두 확률변수 $X$, $Y$가 있을때, 각각의 cdf가 같으면 $X$와 $Y$의 분포가 같다고 볼 수 있음. 

-- 모든 $c$에 대하여 $P(F(X)\leq c) = P(U \leq c )$ 이 성립함을 보이면 된다.

-- 아래와 같이 계산가능 (LHS에서 두번째 등호가 성립하는 이유는 $F$가 순증가함수이기 때문)
-  $RHS=P(U\leq c)=\int_0^c pdf_U(x)dx=\int_0^c 1 dx=c$
-  $LHS=P(F(X)\leq c)=P\big(F^{-1}(F(X))\leq F^{-1}(c)\big)=P(X\leq F^{-1}(c))$

-- 결국에는 $P(X \leq F^{-1}(c)) = c$ 임을 보이면 된다. (아래와 같이 보이면 됨)
- cdf의 정의에 의하여 임의의 $\star$에 대하여 $P(X\leq \star)=F(\star)$ 이 성립. 
-  $\star=F^{-1}(c)$ 를 대입하면 $P(X \leq F^{-1}(c))=F(F^{-1}(c))$ 와 같이 된다. 
- 그런데 $F(F^{-1}(c))=c$ 이므로 $P(X \leq F^{-1}(c))=F(F^{-1}(c))=c$.

"""

# ╔═╡ b61ccc9d-9dd4-486a-b055-90a76fcc0fe7
md"""
**2단계:** $X \overset{d}{=} F^{-1}(U)$ 임을 보이자. 

-- 임의의 $c$ 에 대하여 $P(X\leq c)= P(F^{-1}(U)\leq c)$ 임을 보이면 된다. 
- 임의의 $c$에 대하여 = $c$를 내 마음대로 아무렇게나 잡아도 = 모든 $c$ 에 대하여

-- $RHS = P(F^{-1}(U)\leq c)=P(U \leq F(c))=P(F(X)\leq F(c))=P(X\leq c)$
- 두번째등호: $F$가 순증가함수이므로 성립 
- 세번째등호: $U\overset{d}{=}F(X)$이므로 성립 
- 네번째등호: $F$가 순증가함수이므로 성립 
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
Distributions = "~0.25.107"
Plots = "~1.39.0"
PlutoUI = "~0.7.58"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "bb4c2c3e690fe7dfa615df4a2c4b43a83aa926d6"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "0f748c81756f2e5e6854298f11ad8b2dfae6911a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.0"

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

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

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

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "9c4708e3ed2b799e6124b5673a712dda0b596a9b"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.3.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1fb174f0d48fe7d142e1109a10636bc1d14f5ac2"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.17"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "7c302d7a5fec5214eb8a5a4c466dcf7a51fcf169"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.107"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

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
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

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

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "5b93957f6dcd33fc343044af3d48c215be2562f1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.9.3"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

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

[[deps.Format]]
git-tree-sha1 = "f3cf88025f6d03c194d73f5d13fee9004a108329"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.6"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

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
git-tree-sha1 = "db864f2d91f68a5912937af80327d288ea1f3aee"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

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

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

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
git-tree-sha1 = "3336abae9a713d2210bb57ab484b1e065edd7d23"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.2+0"

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

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cad560042a7cc108f5a4c24ea1431a9221f22c1b"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.2"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e5edc369a598dfde567269dc6add5812cfa10cd5"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.39.3+0"

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
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

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
git-tree-sha1 = "af81a32750ebc831ee28bdaaba6e1067decef51e"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60e3045590bd104a16fefb12836c00c0ef8c7f8c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+0"

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
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

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
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

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
git-tree-sha1 = "862942baf5663da528f66d24996eb6da85218e76"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "ccee59c6e48e6f2edf8a5b64dc817b6729f99eb5"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.39.0"

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
git-tree-sha1 = "71a22244e352aa8c5f0f2adde4150f62368a3f2e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.58"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

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

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

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

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

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
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "cef0472124fab0695b58ca35a77c6fb942fdab8a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.1"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

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
git-tree-sha1 = "3caa21522e7efac1ba21834a03734c57b4611c7e"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.4"
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
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"

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

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "07e470dabc5a6a4254ffebc29a1b3fc01464e105"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.5+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "37195dcb94a5970397ad425b95a9a26d0befce3a"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.0+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

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
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

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
git-tree-sha1 = "1ea2ebe8ffa31f9c324e8c1d6e86b4165b84a024"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+0"

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
# ╟─dbf8cdd6-7815-49f8-9dd6-5987000792ce
# ╟─b7a1c25e-1aeb-4a57-80b5-a9ea4c5c0530
# ╠═20697867-c74a-4485-bb97-c6a31060669a
# ╟─64e29d20-aaf5-4fa0-ad03-b283dac52dce
# ╠═521890de-ab23-11ec-0c2f-2dcaee6dc1bc
# ╠═4e838664-c6ea-4421-99e6-c585fd21cb12
# ╠═61967452-a8e0-46a8-85cb-c0944451c5d2
# ╟─29a2e115-37a7-4f67-a040-a4c1648fdfa6
# ╟─a185e773-f51e-449b-8d40-de819bd4badf
# ╟─78ebbf73-8663-490b-8c0f-58bba2fe2ead
# ╟─36af80d2-bff4-426e-8b44-8f407638941f
# ╟─770f89b0-25ac-42ed-a320-74809e368b44
# ╟─68f0c9b7-25e9-4807-a543-0c7a9a554898
# ╠═51231561-64d1-4545-93d7-980bf77847ae
# ╟─dea01e57-ca18-4a54-9631-8ca3046f4bd8
# ╠═39235f8a-4e2d-4b26-82f4-d95a5c2167ed
# ╟─885941da-8604-4af0-8b69-6d564348f116
# ╟─cb63ae3e-7de4-4eac-aca3-55561864956d
# ╟─2134273b-fca8-4238-b906-f05948fae7d2
# ╠═1da1e41d-73a6-43f4-9d69-3e28b83da8e0
# ╟─caa4ea34-e2ae-455a-88a7-6324ff2cc8ef
# ╟─7620f927-e031-44f5-b22e-bfaa26bdad84
# ╠═eb7572bd-63eb-4792-b982-9f320358a99c
# ╟─1c902aac-d633-4cde-95e1-62e42e0ffe87
# ╟─1e417760-f0c7-4dad-8ffd-003a1ea24f06
# ╟─69324599-c085-45f9-8e6b-359c2d1a163c
# ╠═5d39448b-2e68-475d-a6d4-a06c3b6b347a
# ╟─f3ee3d9c-4163-469a-b96f-8536a413025c
# ╠═fae9c055-f4dd-4d6e-9dc7-ddd6a74fd905
# ╟─aef02ab9-4e12-4030-b428-f4aa36f65ba6
# ╟─6924193e-ea15-42e9-a702-8c6b0a224f35
# ╠═a27efaae-2793-4d8a-8dd6-4900caa80a4a
# ╟─fc1acb74-f856-47c7-b9e1-8d635f0e35ca
# ╟─c2c198a7-4597-49e2-98d7-8776603db7fa
# ╟─cf6a955e-d684-4a16-bc0f-14409f83536a
# ╟─ed777a37-6351-41a1-aed3-86412c9edaac
# ╟─8ab6d3db-1037-4afb-9d30-d16675a1bffd
# ╠═db26ed1c-b01a-4f06-886d-4079cf2c139f
# ╟─7889eeaa-a7b5-4cdd-8b26-c0de2d2981e0
# ╟─0e381090-dac3-48fe-adbb-f25e630d25bf
# ╠═23abf65d-276d-4a5d-ae70-b41255b35b48
# ╟─ac1c2222-04b5-4af0-a55c-723b1ad57dec
# ╟─e4a7a78e-e2c5-4332-ba6e-efe0a87dd4c8
# ╠═9c2f4080-9dd5-4f85-83f7-42bf3c719e6e
# ╟─e234260f-b0bd-4f0e-a5cd-bdc9afad041f
# ╟─cd6507ca-e30c-4c2c-8d4a-11b714db833e
# ╠═b5093fca-2946-41b2-a023-d4abc17c453c
# ╟─d7823cc8-231f-41b4-8a33-175cd42cb2b5
# ╟─35760f05-e335-46cb-baea-83763f3f20d5
# ╟─fbb0395a-9ea9-4eb6-a439-51aa13024a31
# ╟─68dcc6c9-7c52-4d16-b5b0-f0aef055c01a
# ╟─45a26fc2-3c96-4130-83b2-e75c39846264
# ╠═86f22bd5-fb69-46cc-b53e-e1e04e99b5e2
# ╠═f455bad9-640f-4f6d-ab08-04d432681de9
# ╠═6de8f87f-7f99-480a-b5f6-30f0a152f5a1
# ╟─67edeecf-a93a-488d-844b-450866d9b43f
# ╠═5452a29a-2e5a-47c2-804f-7a5cad2d943f
# ╟─72edfa85-7ce4-49ea-a11a-49d714536106
# ╟─b0fefc2a-a6b6-4615-8a94-98a328a09042
# ╟─b66b5d63-9e49-4456-88c3-d47c48870fe9
# ╠═ecbd2a2b-e527-436a-95a1-e62d85cd91b0
# ╟─713adf30-a586-47c7-832c-5b6588203862
# ╟─e88df4b7-5352-49e4-be65-6be16f23f003
# ╟─5d024a81-0244-4348-b6c6-e5e9f70545bd
# ╟─6cf0f4d7-d8f4-4478-8b30-55fbeaeee150
# ╠═bb42b213-ac3d-4156-9d12-445e0702968b
# ╟─413f316d-13a4-4bda-9867-b2439be20f5c
# ╠═1f99cbe2-e472-4916-80bc-ca290ca17f8b
# ╟─1cbdb4a1-0ddb-41e3-8130-ee79005cc193
# ╟─055bcec8-5d53-4888-bf46-50eef3a3e821
# ╠═c565e07f-7bfb-4dfa-9b90-90171d6b5591
# ╟─d31dbd9c-77c6-4c18-802a-c0967ed71cab
# ╟─f23c6456-e364-4658-a71f-6be8388209d0
# ╟─859d98fe-e49e-45d3-9a6f-a6354413ab91
# ╠═b944f877-0984-425b-b35c-09f1d71ffb05
# ╟─6e41df13-d03c-4def-a5b1-2166c945a626
# ╟─f290b98a-31bb-4526-9b56-4543d2108400
# ╟─6bf2ec59-4457-4daf-b143-24b0f8ee4e0e
# ╟─0e963046-9cea-43cc-a4f4-491b59d27845
# ╟─b9a9e616-51d2-409e-863e-2aa9b132a75f
# ╠═de13ff8f-98ec-450c-a609-79b9683afd03
# ╟─cd1a64c8-51f7-423e-8fa2-d67351771482
# ╟─b61ccc9d-9dd4-486a-b055-90a76fcc0fe7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
