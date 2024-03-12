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

# ╔═╡ 2d5df8e4-d187-11ec-1ecb-1d934477ca87
using PlutoUI, Images,LinearAlgebra

# ╔═╡ dbe0d458-a631-4d79-affb-eadf3d9e63dd
md"""
# 5월12일
"""

# ╔═╡ 130b2d27-0a82-4ecc-a27f-cf8631789c14
html"""
<div style="display: flex; justify-content: center;">
<div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
<iframe src=
"
https://www.youtube.com/embed/playlist?list=PLQqh36zP38-zS9I2KLayxNW6jX2dJOhTZ
"
width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ ee7be15f-e730-4137-9822-66e249189882
md"""
## usings
"""

# ╔═╡ af4334ff-940b-4bee-9801-879b6bd69a92
PlutoUI.TableOfContents()

# ╔═╡ 535d7733-5d6f-4e27-a194-b087e87f6850
md"""
## SVD응용: 이미지자료의 압축
"""

# ╔═╡ 198d770c-d0e9-447c-bb42-ef3735ce59d3
md"""
`-` ref: https://math.mit.edu/~gs/linearalgebra/linearalgebra5_7-1.pdf
"""

# ╔═╡ f28cee05-be66-4df4-89f9-4e794ba3f26c
md"""
(예비학습) colorview의 2가지 사용방법
"""

# ╔═╡ 3cca624d-b22b-4a04-aa1f-dbda340b538e
let 
	Mosaic = reshape(rand(3*4*4),(3,4,4))
	#colorview(RGB,Mosaic) # 방법1
	colorview(RGB,Mosaic[1,:,:],Mosaic[2,:,:],Mosaic[3,:,:]) # 방법2
end

# ╔═╡ 59e953f4-40f8-4625-854c-7adc2221164b
md"""
### 프랑스국기
"""

# ╔═╡ 21f8cbae-15a2-41c7-a4c9-6ca77db590d5
md"""
`-` motive: 어떠한 (n,m) 이미지를 표현하기 위해서 사실 항상 nm개의 픽셀이 필요한 것은 아니다. 
"""

# ╔═╡ 47fb9e65-f39a-4668-8c1c-44e0fbe58708
md"""
#### 프랑스 국기의 표현1
"""

# ╔═╡ 5c606ee4-9e75-4b9f-b6d4-2ac5fc3b33f5
France = [
	RGB(0,0,1) RGB(0,0,1) RGB(1,1,1) RGB(1,1,1) RGB(1,0,0) RGB(1,0,0)
	RGB(0,0,1) RGB(0,0,1) RGB(1,1,1) RGB(1,1,1) RGB(1,0,0) RGB(1,0,0)
	RGB(0,0,1) RGB(0,0,1) RGB(1,1,1) RGB(1,1,1) RGB(1,0,0) RGB(1,0,0)
	RGB(0,0,1) RGB(0,0,1) RGB(1,1,1) RGB(1,1,1) RGB(1,0,0) RGB(1,0,0)
]

# ╔═╡ 3a218416-46f4-455d-b2c9-ac284d36cf47
md"""
#### 프랑스 국기의 표현2
"""

# ╔═╡ 01810659-4396-45fb-9856-da477d3d73a4
md"""
`-` France = $\begin{bmatrix} 1\\ 1 \\ 1 \\ 1\end{bmatrix}$ [Blue Blue White White Red Red] 로 표현!
"""

# ╔═╡ f72151ce-595c-482a-8e3d-63ddcb9e42c0
let
	a = [1 ; 1 ; 1 ; 1] # col-vector
	b = [RGB(0,0,1);RGB(0,0,1);RGB(1,1,1);RGB(1,1,1);RGB(1,0,0);RGB(1,0,0)]
	a*b'
end 

# ╔═╡ 2351bf21-5a03-4898-8aad-0569122d5909
md"""
`-` 아래를 전달하지 말고 
```
[
	RGB(0,0,1) RGB(0,0,1) RGB(1,1,1) RGB(1,1,1) RGB(1,0,0) RGB(1,0,0)
	RGB(0,0,1) RGB(0,0,1) RGB(1,1,1) RGB(1,1,1) RGB(1,0,0) RGB(1,0,0)
	RGB(0,0,1) RGB(0,0,1) RGB(1,1,1) RGB(1,1,1) RGB(1,0,0) RGB(1,0,0)
	RGB(0,0,1) RGB(0,0,1) RGB(1,1,1) RGB(1,1,1) RGB(1,0,0) RGB(1,0,0)
]
```
아래를 전달한다면?
```
	[1;1;1;1]*[RGB(0,0,1);RGB(0,0,1);RGB(1,1,1);RGB(1,1,1);RGB(1,0,0);RGB(1,0,0)]' 
```
"""

# ╔═╡ 786d1753-5eea-4268-ad40-d41feeead679
md"""
`-` 그렇다면 $(4*6)*3$ 의 값을 전달하는 대신에 $4+(6*3)$만 전달해도 충분함
"""

# ╔═╡ 5cb92205-62b7-41ba-973d-cadf5d3c1811
4*6*3, 4+6*3 

# ╔═╡ e65a7027-56fb-46a6-aeae-70f5bac5bd07
md"""
#### 프랑스 국기의 표현3: 

`-` 프랑스 국기의 표현2에서 [1;1;1;1]은 모든 RGB에 동일하게 설정하지 않는다면 좀 더 다양한 이미지를 만들 수 있다. 

`-` 즉 Fracne = (${\bf R}$, ${\bf G}$, ${\bf B}$)=(${\boldsymbol a}{\boldsymbol b}^\top$, ${\boldsymbol c}{\boldsymbol d}^\top$, ${\boldsymbol e}{\boldsymbol f}^\top$)의 형태로 표현한다는 의미
"""

# ╔═╡ 85da3736-f645-4961-822b-ffb108ff7f1d
let 
	a = [1,1,1,1] |> float
	b = [0,0,1,1,1,1] |> float
	c = [1,1,0.5,1] |> float
	d = [0,0,1,1,0,0] |> float
	e = [1,1,1,1] |> float
	f = [1,1,1,1,0,0] |> float
	colorview(RGB,a*b',c*d',e*f')
end 

# ╔═╡ 5d9db9f5-d400-472c-99b8-abc0e575530a
md"""
`-` 이렇게 하면 좀더 다양한 이미지를 만들수 있다. 그리고 이전의 프랑스 국기 표현 2에 대비 전달할 픽셀수가 그렇게 많이 증가하지도 않음. 
"""

# ╔═╡ 484a8fbb-bc5b-473f-b5c5-58920377948e
4*6*3, 4+(6*3), (4*3)+(6*3) ## 표현 1,2,3에 사용된 픽셀수 

# ╔═╡ 7185804d-c2f5-42d7-8483-0426410c586c
md"""
### 그리스국기 
"""

# ╔═╡ 99441064-6316-477f-9c9a-a163c127ea25
md"""
`-` 좀 더 복잡한 이미지일 경우
"""

# ╔═╡ cfb98d62-3a0c-45c2-99e2-2c3f144a9c48
begin
	Greece=load(download("https://upload.wikimedia.org/wikipedia/commons/thumb/5/5c/Flag_of_Greece.svg/1200px-Flag_of_Greece.svg.png?20160309091801"))
	Greece=imresize(Greece, (9,13))
end

# ╔═╡ 02e26acb-5370-45eb-9e29-b3ebceb66629
md"""
`-` 프랑스 국기의 경우는 

Fracne = (${\bf R}$, ${\bf G}$, ${\bf B}$)=(${\boldsymbol a}{\boldsymbol b}^\top$, ${\boldsymbol c}{\boldsymbol d}^\top$, ${\boldsymbol e}{\boldsymbol f}^\top$)

와 같이 표현 가능했다. 그리스 국기의 경우는 

Greece = (${\bf R}$, ${\bf G}$, ${\bf B}$)=($\sum_{i=1}^3{\boldsymbol a}_i{\boldsymbol b}_i^\top$, $\sum_{i=1}^3{\boldsymbol c}_i{\boldsymbol d}_i^\top$, $\sum_{i=1}^3{\boldsymbol e}_i{\boldsymbol f}_i^\top$)

로 표현가능하다!!
""" 

# ╔═╡ 6ef4b092-408f-45c8-9ee9-11c1586b467a
md"""
`-` 확인
"""

# ╔═╡ 27e7c978-3902-4a77-872f-8c1faf512cdd
begin 
	a1= [-0.08,-0.42,-0.2,-0.42,-0.08,-0.54,-0.03,-0.54,-0.03]
	a2= [0.07,-0.44,0.72,-0.44,0.07,0.21,0.01,0.21,0.01] 
	a3= [-0.62,-0.15,-0.22,-0.15,-0.62,0.25,0.01,0.25,0.01]
	b1= [-0.2,-0.2,-0.35,-0.2,-0.2,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3,-0.3]
	b2= [0.43,0.43,0.15,0.43,0.43,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17,-0.17]
	b3= [0.15,0.15,-0.92,0.15,0.15,0.09,0.09,0.09,0.09,0.09,0.09,0.09,0.09]
	c1= [-0.21,-0.4,-0.29,-0.4,-0.21,-0.48,-0.18,-0.48,-0.18]
	c2= [0.1,-0.48,0.7,-0.48,0.1,0.13,0.05,0.13,0.05]
	c3= [-0.59,-0.11,-0.13,-0.11,-0.59,0.34,0.13,0.34,0.13]
	d1= -[0.24,0.24,0.34,0.24,0.24,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28]
	d2= [0.41,0.41,0.13,0.41,0.41,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19,-0.19]
	d3= [0.15,0.15,-0.93,0.15,0.15,0.08,0.08,0.08,0.08,0.08,0.08,0.08,0.08]
	e1= [-0.29,-0.36,-0.32,-0.36,-0.29,-0.4,-0.27,-0.4,-0.27]
	e2= [-0.09,0.5,-0.68,0.5,-0.09,-0.08,-0.05,-0.08,-0.05]
	e3= [-0.55,-0.09,-0.11,-0.09,-0.55,0.35,0.24,0.35,0.24]
	f1= -[0.26,0.26,0.31,0.26,0.26,0.28,0.28,0.28,0.28,0.28,0.28,0.28,0.28]
	f2= [-0.4,-0.4,-0.11,-0.4,-0.4,0.21,0.21,0.21,0.21,0.21,0.21,0.21,0.21]
	f3= [0.13,0.13,-0.94,0.13,0.13,0.07,0.07,0.07,0.07,0.07,0.07,0.07,0.07]
	r = 6.84*a1*b1' + 2.56*a2*b2' + 1.38*a3*b3'
	g = 7.80*c1*d1' + 1.75*c2*d2' + 0.82*c3*d3'
	b = 9.36*e1*f1' + 0.88*e2*f2' + 0.37*e3*f3'
	colorview(RGB,r,g,b)
end

# ╔═╡ 95f9f921-6b4d-4ce7-ad23-438468084ed6
md"""
`-` ??? 이걸 어떻게 알았어?! SVD! 
- 이론: 아무리 복잡한 형태의 이미지라도 결국 각 채널(=색깔)은 매트릭스의 형태이므로 $\sum_{i=1}^{n \wedge m}d_iU_iV_i^\top$의 꼴로 표현가능.
- 깨달음: 프랑스 국기의 이미지는 사실 $(size(France)) 이므로 $$\sum_{i=1}^{4}d_iU_iV_i^\top$$로 각 채널이 표현되는데, $$d_2=d_3=d_4=0$$ 이어서 $$d_1U_1V_1^\top$$ 꼴로만 각 채널이 표현된 경우라고 볼 수 있음.
- 정리: 프랑스 국기처럼 단순한 이미지는 각 색깔을 $d_1U_1V_1^\top$로 표현가능하다. 아무리 복잡하더라도 각 색깔을 최소한 $\sum_{i=1}^{n \wedge m}d_iU_iV_i^\top$꼴 로는 표현가능하다. 대부분의 이미지는 프랑스국기처럼 간단한 이미지와 아주 복잡한 이미지의 사이 어딘가쯤에 위치한다. 그리스 국기의 경우 각 채널의 크기가 $(size(Greece)) 인데 $$\sum_{i=1}^{9}d_iU_iV_i^\top$$ 가 아닌 $$\sum_{i=1}^{3}d_iU_iV_i^\top$$ 으로 표현가능하다. 
"""

# ╔═╡ 811dad0d-0695-4efb-bca7-428105d2bd79
md"""
### 하니의 압축 
"""

# ╔═╡ fdc5fa23-7c0a-4ab5-afb5-3769a00b0e77
md"""
(예비학습) 문자열 추가하는 방법, ...의 의미
"""

# ╔═╡ 4ad49abc-759b-46a4-9757-515d02b6215c
push!([1,2,3],4)

# ╔═╡ 60b94fc9-0d39-42de-b2d3-d447bc80495d
push!([1,2,3],[4,5]) # 이것은 에러가 난다.

# ╔═╡ 305a9593-18a2-4762-a68e-08c4088fd47a
push!([1,2,3],4,5) # 이렇게 해야한다.

# ╔═╡ ebdc3a41-c9a8-4752-80a9-e4ae6fb1172a
push!([1,2,3],[4,5]...) # [4,5]의 벡터구조를 그대로 추가하고 싶을 경우 

# ╔═╡ 48faf682-d34e-4df8-a24d-48cdb6d2b7bc
md"""
`-` 의문: 실제로는 각 색깔을 $\sum_{i=1}^{n \wedge m}d_iU_iV_i^\top$로 표현해야 하는 복잡한 이미지인데 만약 $n \wedge m$보다 작은 $k$를 선택하여 $\sum_{i=1}^{k}d_iU_iV_i^\top$로만 표현한다면 어떻게 될까? -> 복잡한 원본과 똑같지는 않겠지만 원본을 간단화 시킨버전의 이미지가 나오지 않을까? 
"""

# ╔═╡ 8245b4f2-4459-4522-8b6a-994c93111652
md"""
`-` 하니를 불러오자. 
"""

# ╔═╡ 2e0107f5-a9ea-4d1d-a8fe-f3db84e59bba
hani = load(download("https://github.com/guebin/SC2022/blob/main/hani.jpeg?raw=true"))

# ╔═╡ 541d4883-1240-44d9-8570-4be52fc1e7ce
md"""
`-` 하니를 텐서로 바꿔보자
"""

# ╔═╡ 23f89b39-3945-413e-881c-4f0c4ce3b7d6
HANI = channelview(hani)

# ╔═╡ 1877f81e-5b05-445e-8597-f2ba265e839a
HANI_SVDs = [HANI[1,:,:], HANI[2,:,:], HANI[3,:,:]] .|> svd

# ╔═╡ 530a5b16-4512-49d3-a967-119024ad2a4a
size(hani) # 3024, 4032 

# ╔═╡ a609e5c5-1eb0-4cc0-9352-29812cb9df73
push!(Array(1:10),(20:10:50)...,100,200,500)

# ╔═╡ 7fd24f6a-3334-421a-85b8-661d0b5e0240
@bind rk Slider(
	push!(Array(1:10),(20:10:50)...,100,200,500),
	default=100,
	show_value=true)

# ╔═╡ 3cf81e78-dfb5-4e3f-bb19-c7c332cdb2e6
function partially_reconstruct(SVDrslt)
	U,d,V = SVDrslt 
	return U[:,1:rk] * Diagonal(d[1:rk]) * V[:,1:rk]'
end

# ╔═╡ 02272643-cafb-43cb-ad4c-45f8674504a6
let 
	HANI_reconstruction =  HANI_SVDs .|> partially_reconstruct
	hani2 = colorview(RGB, HANI_reconstruction...)
	ratio =(3024*rk+4032*rk)*3 / (3024*4032*3)
	md"""
	압축본/원본 = $ratio
	
	[$(hani') $(hani2')]
	"""
end

# ╔═╡ 59b8b39c-086b-43d1-89d2-f3ed971c39af
md"""
- rank=100 정도면 원본과 압축본의 차이를 육안으로 구별하기 힘들다.
"""

# ╔═╡ 853bd9e2-ca6c-4b95-b659-bb8b4d419a85
md"""
`-` SVD를 이용하면 이미지를 압축하여 효율적으로 전달할 수 있다.
- note: SVD는 이미지압축 자체를 위해서 만든 알고리즘이 아니다. (그것에 응용할 뿐)
"""

# ╔═╡ 731681fe-0ac8-4a54-a97a-8100e0c248f7
md"""
## 숙제 
`-` rk=350 일때 압축하니/원본하니의 ratio를 계산하라.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Images = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Images = "~0.25.2"
PlutoUI = "~0.7.38"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

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

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "cf6875678085aed97f52bfc493baaebeb6d40bcb"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.5"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CatIndices]]
deps = ["CustomUnitRanges", "OffsetArrays"]
git-tree-sha1 = "a0f80a09780eed9b1d106a1bf62041c2efc995bc"
uuid = "aafaddc9-749c-510e-ac4f-586e18779b91"
version = "0.2.2"

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

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "75479b7df4167267d75294d14b58244695beb2ac"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.14.2"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "3f1f500312161f1ae067abe07d13b40f78f32e07"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.8"

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

[[deps.ComputationalResources]]
git-tree-sha1 = "52cb3ec90e8a8bea0e62e275ba577ad0f74821f7"
uuid = "ed09eef8-17a6-5b46-8889-db040fac31e3"
version = "0.3.2"

[[deps.CoordinateTransformations]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "681ea870b918e7cff7111da58791d7f718067a19"
uuid = "150eb455-5306-5404-9cee-2592286d6298"
version = "0.6.2"

[[deps.CustomUnitRanges]]
git-tree-sha1 = "1a3f97f907e6dd8983b744d2642651bb162a3f7a"
uuid = "dc8bdbbb-1ca9-579f-8c36-e416f6a65cce"
version = "1.0.2"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

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

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.FFTViews]]
deps = ["CustomUnitRanges", "FFTW"]
git-tree-sha1 = "cbdf14d1e8c7c8aacbe8b19862e0179fd08321c2"
uuid = "4f61f5a4-77b1-5117-aa51-3ab5ef4ef0cd"
version = "0.3.2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "505876577b5481e50d089c1c68899dfb6faebc62"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.6"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "80ced645013a5dbdc52cf70329399c35ce007fae"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.13.0"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78e2c69783c9753a91cdae88a8d432be85a2ab5e"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.0+0"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "1c5a84319923bea76fa145d49e93aa4394c73fc2"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.1"

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "57c021de207e234108a6f1454003120a1bf350c4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.6.0"

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

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageContrastAdjustment]]
deps = ["ImageCore", "ImageTransformations", "Parameters"]
git-tree-sha1 = "0d75cafa80cf22026cea21a8e6cf965295003edc"
uuid = "f332f351-ec65-5f6a-b3d1-319c6670881a"
version = "0.3.10"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

[[deps.ImageDistances]]
deps = ["Distances", "ImageCore", "ImageMorphology", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "7a20463713d239a19cbad3f6991e404aca876bda"
uuid = "51556ac3-7006-55f5-8cb3-34580c88182d"
version = "0.2.15"

[[deps.ImageFiltering]]
deps = ["CatIndices", "ComputationalResources", "DataStructures", "FFTViews", "FFTW", "ImageBase", "ImageCore", "LinearAlgebra", "OffsetArrays", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "TiledIteration"]
git-tree-sha1 = "15bd05c1c0d5dbb32a9a3d7e0ad2d50dd6167189"
uuid = "6a3955dd-da59-5b1f-98d4-e7296123deb5"
version = "0.7.1"

[[deps.ImageIO]]
deps = ["FileIO", "JpegTurbo", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "464bdef044df52e6436f8c018bea2d48c40bb27b"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.1"

[[deps.ImageMagick]]
deps = ["FileIO", "ImageCore", "ImageMagick_jll", "InteractiveUtils", "Libdl", "Pkg", "Random"]
git-tree-sha1 = "5bc1cb62e0c5f1005868358db0692c994c3a13c6"
uuid = "6218d12a-5da1-5696-b52f-db25d2ecc6d1"
version = "1.2.1"

[[deps.ImageMagick_jll]]
deps = ["Artifacts", "Ghostscript_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f025b79883f361fa1bd80ad132773161d231fd9f"
uuid = "c73af94c-d91f-53ed-93a7-00f77d67a9d7"
version = "6.9.12+2"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

[[deps.ImageMorphology]]
deps = ["ImageCore", "LinearAlgebra", "Requires", "TiledIteration"]
git-tree-sha1 = "7668b123ecfd39a6ae3fc31c532b588999bdc166"
uuid = "787d08f9-d448-5407-9aad-5290dd7ab264"
version = "0.3.1"

[[deps.ImageQualityIndexes]]
deps = ["ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "OffsetArrays", "Statistics"]
git-tree-sha1 = "1d2d73b14198d10f7f12bf7f8481fd4b3ff5cd61"
uuid = "2996bd0c-7a13-11e9-2da2-2f5ce47296a9"
version = "0.3.0"

[[deps.ImageSegmentation]]
deps = ["Clustering", "DataStructures", "Distances", "Graphs", "ImageCore", "ImageFiltering", "ImageMorphology", "LinearAlgebra", "MetaGraphs", "RegionTrees", "SimpleWeightedGraphs", "StaticArrays", "Statistics"]
git-tree-sha1 = "36832067ea220818d105d718527d6ed02385bf22"
uuid = "80713f31-8817-5129-9cf8-209ff8fb23e1"
version = "1.7.0"

[[deps.ImageShow]]
deps = ["Base64", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "25f7784b067f699ae4e4cb820465c174f7022972"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.4"

[[deps.ImageTransformations]]
deps = ["AxisAlgorithms", "ColorVectorSpace", "CoordinateTransformations", "ImageBase", "ImageCore", "Interpolations", "OffsetArrays", "Rotations", "StaticArrays"]
git-tree-sha1 = "42fe8de1fe1f80dab37a39d391b6301f7aeaa7b8"
uuid = "02fcd773-0e25-5acc-982a-7f6622650795"
version = "0.9.4"

[[deps.Images]]
deps = ["Base64", "FileIO", "Graphics", "ImageAxes", "ImageBase", "ImageContrastAdjustment", "ImageCore", "ImageDistances", "ImageFiltering", "ImageIO", "ImageMagick", "ImageMetadata", "ImageMorphology", "ImageQualityIndexes", "ImageSegmentation", "ImageShow", "ImageTransformations", "IndirectArrays", "IntegralArrays", "Random", "Reexport", "SparseArrays", "StaticArrays", "Statistics", "StatsBase", "TiledIteration"]
git-tree-sha1 = "03d1301b7ec885b266c0f816f338368c6c0b81bd"
uuid = "916415d5-f1e6-5110-898d-aaa5f9f070e0"
version = "0.25.2"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IntegralArrays]]
deps = ["ColorTypes", "FixedPointNumbers", "IntervalSets"]
git-tree-sha1 = "509075560b9fce23fdb3ccb4cc97935f11a43aa0"
uuid = "1d092043-8f09-5a30-832f-7509e371ab51"
version = "0.1.4"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

[[deps.IntervalSets]]
deps = ["Dates", "Statistics"]
git-tree-sha1 = "eb381d885e30ef859068fce929371a8a5d06a914"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.6.1"

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

[[deps.JLD2]]
deps = ["FileIO", "MacroTools", "Mmap", "OrderedCollections", "Pkg", "Printf", "Reexport", "TranscodingStreams", "UUIDs"]
git-tree-sha1 = "81b9477b49402b47fbe7f7ae0b252077f53e4a08"
uuid = "033835bb-8acc-5ee8-8aae-3f567f8a3819"
version = "0.4.22"

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

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

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

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "c9551dd26e31ab17b86cbd00c2ede019c08758eb"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "76c987446e8d555677f064aaac1145c4c17662f8"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.14"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "e595b205efd49508358f7dc670a940c790204629"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.MetaGraphs]]
deps = ["Graphs", "JLD2", "Random"]
git-tree-sha1 = "2af69ff3c024d13bde52b34a2a7d6887d4e7b438"
uuid = "626554b9-1ddb-594c-aa3c-2596fe9399a5"
version = "0.7.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ded92de95031d4a8c61dfb6ba9adb6f1d8016ddd"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.10"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "e925a64b8585aa9f4e3047b8d2cdc3f0e79fd4e4"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.16"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

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

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.Quaternions]]
deps = ["DualNumbers", "LinearAlgebra", "Random"]
git-tree-sha1 = "b327e4db3f2202a4efafe7569fcbe409106a1f75"
uuid = "94ee1d12-ae83-5a48-8b1c-48b8ff168ae0"
version = "0.5.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RegionTrees]]
deps = ["IterTools", "LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "4618ed0da7a251c7f92e869ae1a19c74a7d2a7f9"
uuid = "dee08c22-ab7f-5625-9660-a9af2021b33f"
version = "0.3.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rotations]]
deps = ["LinearAlgebra", "Quaternions", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "3177100077c68060d63dd71aec209373c3ec339b"
uuid = "6038ab10-8711-5258-84ad-4b1120ba62dc"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays", "Test"]
git-tree-sha1 = "a6f404cc44d3d3b28c793ec0eb59af709d827e4e"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.2.1"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

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

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

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

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

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

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "f90022b44b7bf97952756a6b6737d1a0024a3233"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.5.5"

[[deps.TiledIteration]]
deps = ["OffsetArrays"]
git-tree-sha1 = "5683455224ba92ef59db72d10690690f4a8dc297"
uuid = "06e1c1a7-607b-532d-9fad-de7d9aa2abac"
version = "0.3.1"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─dbe0d458-a631-4d79-affb-eadf3d9e63dd
# ╟─130b2d27-0a82-4ecc-a27f-cf8631789c14
# ╟─ee7be15f-e730-4137-9822-66e249189882
# ╠═2d5df8e4-d187-11ec-1ecb-1d934477ca87
# ╠═af4334ff-940b-4bee-9801-879b6bd69a92
# ╟─535d7733-5d6f-4e27-a194-b087e87f6850
# ╟─198d770c-d0e9-447c-bb42-ef3735ce59d3
# ╟─f28cee05-be66-4df4-89f9-4e794ba3f26c
# ╠═3cca624d-b22b-4a04-aa1f-dbda340b538e
# ╟─59e953f4-40f8-4625-854c-7adc2221164b
# ╟─21f8cbae-15a2-41c7-a4c9-6ca77db590d5
# ╟─47fb9e65-f39a-4668-8c1c-44e0fbe58708
# ╠═5c606ee4-9e75-4b9f-b6d4-2ac5fc3b33f5
# ╟─3a218416-46f4-455d-b2c9-ac284d36cf47
# ╟─01810659-4396-45fb-9856-da477d3d73a4
# ╠═f72151ce-595c-482a-8e3d-63ddcb9e42c0
# ╟─2351bf21-5a03-4898-8aad-0569122d5909
# ╟─786d1753-5eea-4268-ad40-d41feeead679
# ╠═5cb92205-62b7-41ba-973d-cadf5d3c1811
# ╟─e65a7027-56fb-46a6-aeae-70f5bac5bd07
# ╠═85da3736-f645-4961-822b-ffb108ff7f1d
# ╟─5d9db9f5-d400-472c-99b8-abc0e575530a
# ╠═484a8fbb-bc5b-473f-b5c5-58920377948e
# ╟─7185804d-c2f5-42d7-8483-0426410c586c
# ╟─99441064-6316-477f-9c9a-a163c127ea25
# ╠═cfb98d62-3a0c-45c2-99e2-2c3f144a9c48
# ╟─02e26acb-5370-45eb-9e29-b3ebceb66629
# ╠═6ef4b092-408f-45c8-9ee9-11c1586b467a
# ╠═27e7c978-3902-4a77-872f-8c1faf512cdd
# ╟─95f9f921-6b4d-4ce7-ad23-438468084ed6
# ╟─811dad0d-0695-4efb-bca7-428105d2bd79
# ╟─fdc5fa23-7c0a-4ab5-afb5-3769a00b0e77
# ╠═4ad49abc-759b-46a4-9757-515d02b6215c
# ╠═60b94fc9-0d39-42de-b2d3-d447bc80495d
# ╠═305a9593-18a2-4762-a68e-08c4088fd47a
# ╠═ebdc3a41-c9a8-4752-80a9-e4ae6fb1172a
# ╟─48faf682-d34e-4df8-a24d-48cdb6d2b7bc
# ╟─8245b4f2-4459-4522-8b6a-994c93111652
# ╠═2e0107f5-a9ea-4d1d-a8fe-f3db84e59bba
# ╟─541d4883-1240-44d9-8570-4be52fc1e7ce
# ╠═23f89b39-3945-413e-881c-4f0c4ce3b7d6
# ╠═1877f81e-5b05-445e-8597-f2ba265e839a
# ╠═530a5b16-4512-49d3-a967-119024ad2a4a
# ╠═a609e5c5-1eb0-4cc0-9352-29812cb9df73
# ╠═7fd24f6a-3334-421a-85b8-661d0b5e0240
# ╠═3cf81e78-dfb5-4e3f-bb19-c7c332cdb2e6
# ╠═02272643-cafb-43cb-ad4c-45f8674504a6
# ╟─59b8b39c-086b-43d1-89d2-f3ed971c39af
# ╟─853bd9e2-ca6c-4b95-b659-bb8b4d419a85
# ╟─731681fe-0ac8-4a54-a97a-8100e0c248f7
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
