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

# ╔═╡ 28f2f684-2f28-4f44-a41c-70ade81e3927
using Plots, Distributions,PlutoUI,Random

# ╔═╡ 80fe12d2-aa70-4eed-a251-92e481d7e40c
md"""
# 3월29일 
"""

# ╔═╡ 82d1b56f-303d-4412-af8d-dc3638b5186f
html"""
<div style="display: flex; justify-content: center;">
<div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
<iframe src=
"
https://www.youtube.com/embed/playlist?list=PLQqh36zP38-zhn1zsg_GZZtKXLds0caeD
"
width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ 7f9b98a2-a674-4602-8d91-3923e4d7f5e8
md"""
## usings
"""

# ╔═╡ f2990d77-1496-4f7a-996f-a362432f27e6
PlutoUI.TableOfContents()

# ╔═╡ 81fee780-c9e6-42aa-9dca-0637aa6a92f6
Plots.plotly()

# ╔═╡ b9a0a929-debb-4c20-bb35-d664826e8642
md"""
## 지난시간 강의 보충 및 정정내용 
"""

# ╔═╡ ccde90e7-27b4-4c92-8e9a-a19e187535fe
md"""
- 이항분포에서 포아송을 뽑는게 아니고 특정한조건에서 이항분포와 포아송이 비슷해짐 
- `try_until_you_succeed(p)`는 기하분포를 뽑는 함수 맞아요! 
- `sample size`를 의미하는 기호를 `N`으로 통일. (`n`은 이항분포의 $n$)
- 포아송의 역수로 지수분포를 못 뽑는 이유? (1) 0이 나오니까 역수를 취할수 없음 (2) 지수분포는 임의의 양의값이 가능해야 하는데 $\frac{1}{0.5}, \frac{1}{1.5}$와 같은 값은 포아송의 역수로 구할 수 없음. 
- 강의노트에 제목수준 및 내용을 조금씩 수정함.  
"""

# ╔═╡ 66f787b4-5cda-4883-be6c-b3e9ba5ea3b7
md"""
## 지수분포
"""

# ╔═╡ 31a07e09-fac8-4082-be3c-9d83fd0cb09d
md"""
### How to generate it 
"""

# ╔═╡ 8546e329-714f-4a5f-a7e4-f2dfbf55496d
md"""
(방법1) 모듈 
"""

# ╔═╡ 6ea24d60-1ddc-4462-af52-5009136051f8
md"""
(방법2) 포아송프로세스를 이용
"""

# ╔═╡ 2f7acc36-892b-4368-ac9b-3ae42e20a681
md"""
(방법3) inverse cdf method
"""

# ╔═╡ 662f2f6b-7c45-4a3a-b936-a0cf95a4c908
md"""
##### Inverse cdf method motive 
"""

# ╔═╡ 4f98ca0f-580e-46b1-a279-0bfda3d40dd2
md"""
`-` 아래와 같은 2개의 지수분포를 고려해보자. 
- 평균이 1인 지수분포 
- 평균이 5인 지수분포 
"""

# ╔═╡ 53fc9117-db8c-424f-9098-e2b887240a72
md"""
`-` 두 지수분포에 대한 cdf는 아래와 같다. 
"""

# ╔═╡ be2dee73-c901-4f3f-aa96-c62a886b7167
let
	p1 = plot(x -> -exp(-x)+1,0,20, title="평균=1")
	p2 = plot(x -> -exp(-x/5)+1,0,20, title="평균=5")
	plot(p1, p2, layout=(1,2))
end

# ╔═╡ 24a3f936-1cda-4d4b-bde7-1dbe81e44909
md"""
`-` cdf의 y축에서 랜덤변수를 발생시킨뒤에 $\rightarrow \downarrow$와 같이 이동하여 x축에 내린다고 생각해보자. 
- 왼쪽: 대부분 5이하에 떨어짐
- 오른쪽: 대략 60퍼정도만 5이하에 떨어짐 
"""

# ╔═╡ dbe9914c-da72-41e5-895f-24394097bf22
md"""
`-` 구현 
"""

# ╔═╡ ff22c6d4-7808-4452-8597-f307dd099e7e
let
	Finv(x) = -log(1-x) # 평균이 1인 지수분포 cdf의 역함수 
	Ginv(x) = -5log(1-x) # 평균이 5인 지수분포 cdf의 역함수
	u = rand(5) # 5개의 샘플을 유니폼에서 추출
	p1 = plot(x -> -exp(-x)+1,0,20, title="평균=1")
	scatter!([0,0,0,0,0],u) 
	scatter!(Finv.(u),[0,0,0,0,0])
	p2 = plot(x -> -exp(-x/5)+1,0,20, title="평균=5")
	scatter!([0,0,0,0,0],u)
	scatter!(Ginv.(u),[0,0,0,0,0])
	plot(p1, p2, layout=(1,2))	
end

# ╔═╡ 8507fb19-40e9-4bf7-95f5-c107246723e2
md"""
- 빨간색: 균등분포 
- 초록색: 이게 지수분포 같은데? 
"""

# ╔═╡ e6db3918-4198-4a46-aa5b-663d6c959369
md"""
##### inverse cdf method 알고리즘 정리 

확률변수 $X_1,X_2,\dots,X_n \overset{iid}{\sim} F$ 를 생성하고 싶다면? 

1.  균등분포에서 $n$개의 난수를 독립적으로 생성한다. 이를 $U_1,U_2, \dots U_n$ 이라고 하자. 
2.  $X_1 = F^{-1}(U_1),X_2=F^{-1}(U_2),\dots, X_n=F^{-1}(U_n)$ 이라고 놓는다. 
"""

# ╔═╡ 1480a122-73c7-4704-9299-d21fd18feedf
md"""
`-` 예제1: inverse cdf 를 이용하여 평균이 1인 지수분포 10000개를 생성하여 보자. 
"""

# ╔═╡ 5a95784f-89a8-42cc-95a5-9b0310a23ca7
md"""
(풀이)
"""

# ╔═╡ f53186b7-7fbc-4fe7-8f65-223e681e87e5
rand(10000)

# ╔═╡ f68bd58a-8b00-48ad-8ce1-ba2d19034b95
let
	p1= histogram(rand(10000) .|> x -> -log(1-x))
	p2= histogram(rand(Exponential(1),10000))
	plot(p1,p2,layout=(2,1))  
end

# ╔═╡ bc64fe93-e916-49ce-8a1a-084098c2323f
md"""
### 지수분포의 무기억성 
"""

# ╔═╡ 595862ae-0ee5-47fc-b388-c3986d2c5357
md"""
`-` 이론: $X \sim Exp(\lambda) \Rightarrow$ 모든 $t,s > 0$ 에 대하여 $P(X>t)=P(X>t+s|X>s)$가 성립
"""

# ╔═╡ 1d8dfc66-6b9d-4d8f-b06e-719897cdba3a
md"""
`-` 개념: 

- 이해를 위해서 $t=1,s=9$ 대입 => $P(X>1)=P(X>10 | X>9)$
- 좌변: 시간을 1 기다려서 이벤트가 발생안할 확률 
- 우변: 시간을 9 기다렸는데 이벤트가 발생안했음 -> 시간 10을 기다려서 이벤트가 발생안할 확률 
- 예를들어서 $\lambda=0.1$ 이라면 한번 이벤트 발생하는데 평균 시간10이 걸린다는 의미임. => (1) 좌변은 이제 시간1기다림 (2) 우변은 시간 9를 기다림. 곧 "약속된" 시간 10이 완성됨 => 우변이 더 확률이 크지 않을까? => 아니라는거에요!
"""

# ╔═╡ e9fa46b8-1e1b-4cdf-99a8-77b0dc18a29b
md"""
`–` 이해: 지수분포의 근본? 포아송 프로세스 
- 엄청 짧은 시간 
- 엄청 작은 확률 
- 엄청 많은 베르누이 시행이 "독립적"으로 수행 -> 지금까지 실패했다고 해서 이후에 성공확률이 높아지는건 아님. 
- 우변: 이미 시간9동안 무수히 많은 독립적인 베르누이 시행을 놓친상태임. 그 이후의 시행은 모두 독립이므로 좌변의 확률보다 더 크다고 볼 수 없음. 
"""

# ╔═╡ b5d6e46f-4884-4d7d-a626-1e357fb5716b
md"""
`-` 시뮬레이션으로 확인해보자. 
"""

# ╔═╡ bfa27b33-3a88-4ce4-a43f-3a5a7a3fdb88
md"t= $@bind t Slider(0.01:0.01:5,show_value=true)"

# ╔═╡ 235f35dd-e81f-4aff-a0ae-0355f015d7d1
md"s= $@bind s Slider(0.01:0.01:5,show_value=true)"

# ╔═╡ 91dd34f0-7710-4768-93fd-d50d2ca8dd2e
let
	N=500000
	X = rand(Exponential(1),N) 
	lhs= length(X[X.>t])/N
	rhs= length(X[X.>t+s])/length(X[X.>s])
	md"""
	-  $t=$ $t
	-  $s=$ $s
	-  $P(X>t)=$ $(lhs)
	-  $P(X>t+s|X>s)=$ $(rhs)
	"""
end

# ╔═╡ 17a38ef5-7fa0-4471-9da0-bd910f594b48
md"""
`–` 무기억성 = 과거는 중요하지 않음! => $P(X>1)=P(X>2|X>1)=P(X>3|X>2)=...$ 
"""

# ╔═╡ 9f8e6b38-edf3-428d-b065-bd4f9ea6564e
let 
	N = 5000
	X = rand(Exponential(1),N)
	prob1 = length(X[X.>1])/N
	prob2 = length(X[X.>2])/length(X[X.>1])
	prob3 = length(X[X.>3])/length(X[X.>2])
	prob4 = length(X[X.>4])/length(X[X.>3])
	prob5 = length(X[X.>5])/length(X[X.>4])
	md"""
	-  $P(X>1)$= $(prob1)
	-  $P(X>2|X>1)$= $(prob2) 
	-  $P(X>3|X>2)$= $(prob3)
	-  $P(X>4|X>3)$= $(prob4)
	-  $P(X>5|X>4)$= $(prob5)
	"""
end

# ╔═╡ d2aeff2e-9a70-4a53-a36d-fb9e8d19fdab
md"""
### 몬테카를로 적분 
"""

# ╔═╡ 6eb0b3f4-369e-4197-b8aa-13d745d76078
md"""
##### 예제1: 아래를 계산하라. 
$\int_0^{\infty} x e^{-x}dx =?$
"""

# ╔═╡ 7ef592bd-76f7-49c2-ad05-1a3bac201b00
md"""
(손풀이)
$\int_0^\infty x e^{-x}dx =??=1$
"""

# ╔═╡ b015708a-b199-4999-9f52-088baf8eff8e
md"""
(손풀이2)
$\int_0^\infty x e^{-x}dx =\int_0^\infty x \times e^{-x}dx$은 $\lambda=1$인 지수분포의 평균이다. 따라서 답은 1.
"""

# ╔═╡ 798d8f8e-4447-4f1f-a394-1e9819796853
md"""
(컴퓨터를 이용한 풀이)
"""

# ╔═╡ 7fae9ab3-86d1-4400-ab12-e50d4fdcb480
rand(Exponential(1),10000) |> mean

# ╔═╡ d44bf5cd-7b8e-4e18-86fb-25e0e413e86e
md"""
##### 예제2: 아래를 계산하라. 
$\int_0^{\infty} x^2 e^{-x}dx=?$
"""

# ╔═╡ 3be1e916-fcea-4f09-9b33-911b1e24afb2
md"""
(컴퓨터를 이용한 풀이)
"""

# ╔═╡ fdd7221c-6550-4375-ba81-2678769b2c1f
rand(Exponential(1),10000) .|> (x->x^2) |> mean

# ╔═╡ 1858c174-a3a0-4a71-883a-265e571b05b3
md"""
- 분산 = 제곱의평균 - 평균의제곱 => 제곱의평균 = 분산 + 평균의제곱 = 1+1²
"""

# ╔═╡ 1117b745-9711-4a7f-8569-51def002998a
md"""
##### 예제3: 아래를 계산하라. 
$\int_0^1 e^{-x} dx=?$
"""

# ╔═╡ 3fec3dba-b8ab-4186-acca-641e593f00ec
md"""
(컴퓨터를 이용한 풀이)
"""

# ╔═╡ fc50de2c-82d7-49d7-9a2b-f3a16fabd05e
md"""
구하는것은 $\int_0^{\infty} I(x \leq 1)e^{-x}dx=E(I(X\leq 1)),~ X\sim Exp(1)$ 로 해석 가능하므로 
"""

# ╔═╡ 69bb880f-4a28-430c-bb99-2668f2962b59
md"""
- 단, $I(x \in A )= \begin{cases} 1 & x\in A \\ 0 & o.w \end{cases}$
"""

# ╔═╡ b8f80675-dacf-44c6-9f21-43a90bdf14c7
let
	X = rand(Exponential(1),10000)
	f = x -> (x<=1)
	X .|> f |> mean 
end

# ╔═╡ ae579f27-6f9e-4998-97d3-9db8ab01d7e4
md"""
(컴퓨터를 이용한 풀이2)
"""

# ╔═╡ 808fecdb-9768-45b7-9c20-2744e33ed8e7
md"""
구하는것은 $E(I(X\leq 1))=P(X\leq 1),~ X\sim Exp(1)$ 로 해석 가능하므로 
"""

# ╔═╡ 35aa02e1-6556-4c79-8ee8-10f1c08eddec
let
	N = 10000
	X = rand(Exponential(1),N)
	length(X[X.<=1])/N
end

# ╔═╡ a5d94a21-12b1-47b7-935d-0d937de14550
md"""
(컴퓨터를 이용한 풀이3)
"""

# ╔═╡ eb6ceb7d-f64d-4c6e-a31e-df73ee9f33c3
md"""
구하는것은 $\int_0^1 e^{-x} dx=\int_0^1 e^{-x}\times 1 dx=E(e^{-X}), ~ X\sim U(0,1)$ 로 해석 가능하므로 
"""

# ╔═╡ d2131950-e408-45a7-80a0-002d17d505ce
let 
	N = 10000
	X = rand(N)
	exp.(-X) |> mean
end 

# ╔═╡ 6cfdd6e2-8f31-4a3e-8f67-0651480310b5
md"""
`-` 참고로 세번째 풀이를 응용하면 아래와 같이 적분구간이 유한한 경우에는 쉽게 그 값을 계산할 수 있다. 
"""

# ╔═╡ ebfb2011-0c95-49d8-81f0-fffe45318c18
md"""
$\int_0^{\pi}\sin(x)dx=\int_0^{\pi}\sin(x)\pi \times \frac{1}{\pi}dx=E(\sin(X)\pi),~ X \sim U(0,\pi)$
"""

# ╔═╡ 80fbc0bd-a137-4883-afb7-cfc44954e33e
rand(10000)*π # U(0,2π)

# ╔═╡ 5f5b0fe7-cfb2-4190-9ac1-d76cca5f77ff
mean(rand(10000)*π .|> x-> sin(x)*π)

# ╔═╡ ceee1a6f-f2ec-4af8-a114-29bc58e59d8e
-cos(π)+cos(0)

# ╔═╡ 5d0c9bff-a6c7-4c76-8b5e-71bb3197ba50
md"""
### 박스뮬러변환
"""

# ╔═╡ 49561af0-ac59-4fc2-ac52-13e4239c0fe0
md"""
`-` 이론: $\begin{cases} R^2/2 \sim Exp(1) \\ \Theta \sim U(0,2\pi) \end{cases} \Rightarrow \begin{bmatrix} R\cos \Theta \\ R \sin \Theta \end{bmatrix} \sim N\left (\begin{bmatrix} 0 \\ 0 \end{bmatrix}, \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix}\right)$ 
"""

# ╔═╡ 42ec1e91-1152-47b2-a43a-55a4f7b19743
md"""
- 의미?
"""

# ╔═╡ ba38910a-8b0e-405b-a909-2ba77ca0ba13
@bind i Slider(1:5000)

# ╔═╡ 809c3606-8c0c-400c-a698-5d905bf0ebaa
let 
	Random.seed!(43052)
	X=randn(5000)
	Y=randn(5000)
	plot(X,Y,seriestype=:scatter,alpha=0.1)
	plot!([0,X[i]],[0,Y[i]])
end

# ╔═╡ 26512f21-7de0-4e50-9a76-bf84e64492a4
let 
	N=5000
	X = rand(Normal(0,1),N)
	Y = rand(Normal(0,1),N)
	p1= (X.^2 + Y.^2)./2 |> histogram
	p2= rand(Exponential(1),N) |> histogram
	plot(p1,p2,layout=(2,1))
end 

# ╔═╡ e2b0ec64-0f07-4aeb-a7bb-dcfe2d20b300
md"""
`-` 그렇다면 이것을 응용하여 정규분포를 생성할 수 있겠다. 
"""

# ╔═╡ b4c8ba40-9409-4ebb-b63c-cae56eebdf26
let 
	N= 10000000
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

# ╔═╡ a62728c4-e91a-4345-b0cb-ce1e243f0c81
md"""
`-` inverse cdf 기법과 합치면 아래와 같이 정리가능하다. 

$\begin{cases}
X=\sqrt{-2\log(1-U_1)} \cos(2\pi U_2) \\ 
Y=\sqrt{-2\log(1-U_1)} \sin(2\pi U_2) 
\end{cases},~ U_1,U_2 \overset{iid}{\sim} U(0,1)$
"""

# ╔═╡ 2a7d5045-5391-4769-9d94-3e5deb92924a
let 
	N=10000000
	U₁ = rand(N)
	U₂ = rand(N)
	X = @. √(-2log(1-U₁))*cos(2π*U₂)
	p1 = histogram(X) 
	p2 = histogram(randn(N))
	plot(p1,p2,layout=(2,1),xlim=(-5,5))
end 

# ╔═╡ 0fbbfaea-411d-4cf7-82b4-113038fab29b
md"""
(보충학습) 브로드캐스팅 매크로 
"""

# ╔═╡ 89595f20-6287-4637-a832-b3a2a438b499
[1,2] .+ 1 .+ 1 .+ 1 .+ 1 .+ 1

# ╔═╡ 71114a46-b1c5-4cb1-8ce1-c69200892eb4
@. [1,2] + 1 + 1 + 1 + 1 + 1

# ╔═╡ 37265165-a917-4901-89fc-da994967bbed
md"""
### $\lambda$에 따른 포아송과 지수분포의 히스토그램변화 관찰 
"""

# ╔═╡ 959b6759-91b6-41b9-a67c-43a09eab0857
lambda = @bind λ Slider(0.1:0.01:50,default=10)

# ╔═╡ 1fa64776-6859-4698-bca8-4286bdc01128
let 
	p1=histogram(rand(Poisson(λ),10000),title="포아송, λ=$λ")
	p2=histogram(rand(Exponential(1/λ),10000),title="지수분포, λ=$λ")
	plot(p1,p2,layout=(2,1))
end 

# ╔═╡ 57da91d7-17e5-4a22-b493-62f9d590151d
md"""
### 숙제 

(1) inverse cdf를 이용하여 평균이 1인 지수분포를 발생 

(2) 서로독립인 2개의 정규분포를 이용하여 평균이 1인 지수분포를 발생 

(3) 위의 2개의 히스토그램을 비교 
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
Distributions = "~0.25.52"
Plots = "~1.27.3"
PlutoUI = "~0.7.37"
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

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

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
git-tree-sha1 = "96b0bc6c52df76506efc8a441c6cf1adcb1babc4"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.42.0"

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
git-tree-sha1 = "c43e992f186abaf9965cc45e372f4693b7754b22"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.52"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ae13fcbc7ab8f16b0856729b050ef0c446aa3492"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.4+0"

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
git-tree-sha1 = "9f836fb62492f4b0f0d3b06f55983f2704ed0883"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.64.0"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a6c850d77ad5118ad3be4bd188919ce97fffac47"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.64.0+0"

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

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "65e4589030ef3c44d3b90bdc5aac462b4bb05567"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.8"

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
git-tree-sha1 = "4f00cc36fede3c04b8acf9b2e2763decfdcecfa6"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.13"

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
git-tree-sha1 = "58f25e56b706f95125dcb796f39e1fb01d913a71"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.10"

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
git-tree-sha1 = "85b5da0fa43588c75bb1ff986493443f821c70b7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.2.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "5f6e1309595e95db24342e56cd4dabd2159e0b79"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.27.3"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "bf0a1121af131d9974241ba53f601211e9303a9e"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.37"

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
git-tree-sha1 = "c3d8ba7f3fa0625b062b82853a7d5229cb728b6b"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.2.1"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "25405d7016a47cf2bd6cd91e66f4de437fd54a07"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.16"

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
# ╟─80fe12d2-aa70-4eed-a251-92e481d7e40c
# ╟─82d1b56f-303d-4412-af8d-dc3638b5186f
# ╟─7f9b98a2-a674-4602-8d91-3923e4d7f5e8
# ╠═28f2f684-2f28-4f44-a41c-70ade81e3927
# ╠═f2990d77-1496-4f7a-996f-a362432f27e6
# ╠═81fee780-c9e6-42aa-9dca-0637aa6a92f6
# ╟─b9a0a929-debb-4c20-bb35-d664826e8642
# ╟─ccde90e7-27b4-4c92-8e9a-a19e187535fe
# ╟─66f787b4-5cda-4883-be6c-b3e9ba5ea3b7
# ╟─31a07e09-fac8-4082-be3c-9d83fd0cb09d
# ╟─8546e329-714f-4a5f-a7e4-f2dfbf55496d
# ╟─6ea24d60-1ddc-4462-af52-5009136051f8
# ╟─2f7acc36-892b-4368-ac9b-3ae42e20a681
# ╟─662f2f6b-7c45-4a3a-b936-a0cf95a4c908
# ╟─4f98ca0f-580e-46b1-a279-0bfda3d40dd2
# ╟─53fc9117-db8c-424f-9098-e2b887240a72
# ╠═be2dee73-c901-4f3f-aa96-c62a886b7167
# ╟─24a3f936-1cda-4d4b-bde7-1dbe81e44909
# ╟─dbe9914c-da72-41e5-895f-24394097bf22
# ╠═ff22c6d4-7808-4452-8597-f307dd099e7e
# ╟─8507fb19-40e9-4bf7-95f5-c107246723e2
# ╟─e6db3918-4198-4a46-aa5b-663d6c959369
# ╟─1480a122-73c7-4704-9299-d21fd18feedf
# ╟─5a95784f-89a8-42cc-95a5-9b0310a23ca7
# ╠═f53186b7-7fbc-4fe7-8f65-223e681e87e5
# ╠═f68bd58a-8b00-48ad-8ce1-ba2d19034b95
# ╟─bc64fe93-e916-49ce-8a1a-084098c2323f
# ╟─595862ae-0ee5-47fc-b388-c3986d2c5357
# ╟─1d8dfc66-6b9d-4d8f-b06e-719897cdba3a
# ╟─e9fa46b8-1e1b-4cdf-99a8-77b0dc18a29b
# ╟─b5d6e46f-4884-4d7d-a626-1e357fb5716b
# ╠═bfa27b33-3a88-4ce4-a43f-3a5a7a3fdb88
# ╠═235f35dd-e81f-4aff-a0ae-0355f015d7d1
# ╟─91dd34f0-7710-4768-93fd-d50d2ca8dd2e
# ╟─17a38ef5-7fa0-4471-9da0-bd910f594b48
# ╠═9f8e6b38-edf3-428d-b065-bd4f9ea6564e
# ╟─d2aeff2e-9a70-4a53-a36d-fb9e8d19fdab
# ╟─6eb0b3f4-369e-4197-b8aa-13d745d76078
# ╟─7ef592bd-76f7-49c2-ad05-1a3bac201b00
# ╟─b015708a-b199-4999-9f52-088baf8eff8e
# ╟─798d8f8e-4447-4f1f-a394-1e9819796853
# ╠═7fae9ab3-86d1-4400-ab12-e50d4fdcb480
# ╟─d44bf5cd-7b8e-4e18-86fb-25e0e413e86e
# ╟─3be1e916-fcea-4f09-9b33-911b1e24afb2
# ╠═fdd7221c-6550-4375-ba81-2678769b2c1f
# ╟─1858c174-a3a0-4a71-883a-265e571b05b3
# ╟─1117b745-9711-4a7f-8569-51def002998a
# ╟─3fec3dba-b8ab-4186-acca-641e593f00ec
# ╟─fc50de2c-82d7-49d7-9a2b-f3a16fabd05e
# ╟─69bb880f-4a28-430c-bb99-2668f2962b59
# ╠═b8f80675-dacf-44c6-9f21-43a90bdf14c7
# ╟─ae579f27-6f9e-4998-97d3-9db8ab01d7e4
# ╟─808fecdb-9768-45b7-9c20-2744e33ed8e7
# ╠═35aa02e1-6556-4c79-8ee8-10f1c08eddec
# ╟─a5d94a21-12b1-47b7-935d-0d937de14550
# ╟─eb6ceb7d-f64d-4c6e-a31e-df73ee9f33c3
# ╠═d2131950-e408-45a7-80a0-002d17d505ce
# ╟─6cfdd6e2-8f31-4a3e-8f67-0651480310b5
# ╟─ebfb2011-0c95-49d8-81f0-fffe45318c18
# ╠═80fbc0bd-a137-4883-afb7-cfc44954e33e
# ╠═5f5b0fe7-cfb2-4190-9ac1-d76cca5f77ff
# ╠═ceee1a6f-f2ec-4af8-a114-29bc58e59d8e
# ╟─5d0c9bff-a6c7-4c76-8b5e-71bb3197ba50
# ╟─49561af0-ac59-4fc2-ac52-13e4239c0fe0
# ╠═42ec1e91-1152-47b2-a43a-55a4f7b19743
# ╠═ba38910a-8b0e-405b-a909-2ba77ca0ba13
# ╠═809c3606-8c0c-400c-a698-5d905bf0ebaa
# ╠═26512f21-7de0-4e50-9a76-bf84e64492a4
# ╟─e2b0ec64-0f07-4aeb-a7bb-dcfe2d20b300
# ╠═b4c8ba40-9409-4ebb-b63c-cae56eebdf26
# ╟─a62728c4-e91a-4345-b0cb-ce1e243f0c81
# ╠═2a7d5045-5391-4769-9d94-3e5deb92924a
# ╟─0fbbfaea-411d-4cf7-82b4-113038fab29b
# ╠═89595f20-6287-4637-a832-b3a2a438b499
# ╠═71114a46-b1c5-4cb1-8ce1-c69200892eb4
# ╟─37265165-a917-4901-89fc-da994967bbed
# ╠═959b6759-91b6-41b9-a67c-43a09eab0857
# ╠═1fa64776-6859-4698-bca8-4286bdc01128
# ╟─57da91d7-17e5-4a22-b493-62f9d590151d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
