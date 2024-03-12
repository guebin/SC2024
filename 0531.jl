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

# ╔═╡ b868981c-dffb-4165-8e69-97889bc25312
using LinearAlgebra,PlutoUI

# ╔═╡ 27d01165-2aa4-4b86-93ba-c7109636e692
md"""
# 5월31일
"""

# ╔═╡ 1d6ec17b-eca2-4289-8cfe-04a67bcac2c6
html"""
<div style="display: flex; justify-content: center;">
<div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
<iframe src=
"
https://www.youtube.com/embed/playlist?list=PLQqh36zP38-zj2V6OUIAP5NmybY2LUszL
"
width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ d3ccc0c6-e092-11ec-3296-b5a32d51eade
md"""
## usings
"""

# ╔═╡ 26300fa0-9f73-4874-9e60-d7b82f05cc71
PlutoUI.TableOfContents()

# ╔═╡ 12aa1c0d-8f14-48cf-b244-c3fac884d676
md"""
## Spectral theorem
"""

# ╔═╡ 43a2c50f-0427-4c2c-af3e-c83d92c8e82c
md"""
`-` 이론: 모든 원소가 실수인 대칭행렬 (1) 대각화가능하고 (2) 실수의 고유값을 가지며 (3) 고유벡터행렬이 직교행렬이다. 즉 모든 원소가 실대칭행렬 ${\bf A}$는 아래와 같이 표현할 수 있다. 

${\bf A} = {\bf \Psi}{\bf \Lambda}{\bf \Psi}^\top$

단 ${\bf \Lambda}$의 모든 대각선 원소는 실수이다. 

- note: 이 정리를 스펙트럼정리 (spectral theorem) 혹은 주축정리 (principal axis theorem) 이라고 한다. 
- note: 이러한 매트릭스 ${\bf A}$를 직교대각화가능 (orthogonally diagonalizable) 이라고 부른다. 
"""

# ╔═╡ 742c33c4-3bd5-4388-a9dd-4f9341eeb3e4
md"""
(예제1) 0행렬은 모든 원소가 실수이며 대칭행렬이다. 따라서 직교대각화 가능하다. 
"""

# ╔═╡ 9b52b105-40b4-4b77-b274-2d0a6914985a
let 
	A = [0 0 ; 0 0]
	eigen(A)
end 

# ╔═╡ f1c60bf5-2052-4bf3-8e3b-ddb0cbd14b39
md"""
(예제2) 아무렇게나 실수인 대칭행렬을 만들기만하면 직교대각화가 가능하다.
"""

# ╔═╡ 0b882110-64b9-49c7-8ffd-66b95f0b8a5d
let 
	A = [-33 1 
		 1 -0.22]
	λ,Ψ = eigen(A)
	# Ψ*Ψ'
end 

# ╔═╡ bd476cf8-9946-4911-bee9-5170aabd9ca7
md"""
(예제3) 대칭행렬이지만 모든 원소가 실수가 아닌 경우에는 대각화 불가능 할 수도 있다. 
"""

# ╔═╡ 2430ceb6-1154-4345-8487-05bd9a8136e2
let 
	A = [2im 1 
		1 0]
	λ,Ψ = eigen(A)
	#rank(Ψ)
end 

# ╔═╡ 05af5915-30a6-49fd-850a-352a41036984
md"""
- 이때 $\lambda=i$ 이고 $\psi_1=\frac{1}{\sqrt{2}} \begin{bmatrix} 1 \\ -i \end{bmatrix}$, $\psi_2 = \frac{1}{\sqrt{2}} \begin{bmatrix} i \\ 1 \end{bmatrix}$ 이다. 
-  $\psi_1 i = \psi_2$임을 관찰하라. (따라서 고유벡터행렬이 full rank가 아니다)
- 즉 이 행렬은 대각화가 불가능한 행렬이다. 
"""

# ╔═╡ 9022b592-6ab0-40de-b96b-aecb56c663a5
let 
	ψ1= 1/√2 .* [1,-im]
	ψ2= 1/√2 .* [im,1]
	Ψ = [ψ1 ψ2]
	rank(Ψ)
	#inv(Ψ)
end 

# ╔═╡ f080d7b0-4c1f-411e-9e86-d422dd77d2e7
md"""
## Positive definite matrix
"""

# ╔═╡ c226e8e1-e97c-4039-8baa-f82a00254917
md"""
### Positive definite matrix 정의 
"""

# ╔═╡ 031acc24-0027-4db0-9270-b63a026e9575
md"""
`-` 정의: 어떠한 매트릭스 ${\bf A}_{n\times n}$가 모든 non-zero vector ${\bf y}_{n \times 1}$에 대하여 

$${\bf y}^\top {\bf A}{\bf y} > 0$$

을 만족하면 ${\bf A}$를 positive definite matrix라고 부른다. 

- note: 참고로 등호가 있는 경우는 positive semidefinite matrix라고 부른다. 
"""

# ╔═╡ eac996bd-5458-40eb-a82d-88f3aac35033
md"""
(예제1) 예를들면 ${\bf A}=\begin{bmatrix} 2 & 0 \\ 2 & 2 \end{bmatrix}$는 positive definite 의 정의를 만족한다. 
"""

# ╔═╡ 4d737613-177c-45e8-abe7-597ffd863049
@bind y1 Slider(-50:0.1:50, show_value=true)

# ╔═╡ df004a39-17f0-47eb-9725-d46dd844f205
@bind y2 Slider(-50:0.1:50, show_value=true)

# ╔═╡ 4eabb492-2f64-4eb1-9041-a8d29eb40ff8
let 
	A = [2 0; 2 2]
	y=[y1,y2]
	y'*A*y
end

# ╔═╡ 5b94d293-7533-4f4d-8c71-0a2c19e94117
md"""
(why?) 

${\bf y}'{\bf A}{\bf y}=(y_1+y_2)^2+y_1^2+y_2^2>0$
"""

# ╔═╡ 1bcbd243-7088-4133-a8a9-6dab8d37001c
md"""
- note: 참고로 이때 ${\bf A}$는 대칭행렬이 아니다. 
"""

# ╔═╡ 94b5adf6-cbe5-4052-a2bd-dfabac6d7406
md"""
(예제2) ${\bf A}=\begin{bmatrix} 1 & 1 \\ -1 & 1 \end{bmatrix}$ 는 positive definite의 정의를 만족한다. 
"""

# ╔═╡ 6cc96118-41ac-4807-8f88-4ae86e491c79
let
	A = [1 1 
		-1 1]
	y=[y1,y2]
	y'*A*y
	#eigen(A)
end 

# ╔═╡ dc456627-72c5-4acf-89af-aa1ba01d2ba4
md"""
- note: 참고로 이때 ${\bf A}$의 고유값이 양수인 것은 아니다. 
"""

# ╔═╡ fb480817-66f5-468b-b6c6-9bf145fc86c4
md"""
### Positive definite and symmetric matrix
"""

# ╔═╡ 3e069b39-28fb-4404-b292-8a54c2a25a8a
md"""
`-` 관찰1: 어떠한 실대칭행렬 ${\bf A}_{n\times n}$의 모든 고유값이 양수이면 ${\bf A}$는 positive definite matrix가 된다. (즉 모든 non-zero vector ${\bf y}$에 대해서 ${\bf y}^\top {\bf A}{\bf y}>0$을 만족한다.)
"""

# ╔═╡ 7befd7ef-7ab9-4a00-a9d8-52cde80c23e6
md"""
(왜?) 
"""

# ╔═╡ 71d85e9f-a649-4584-8cbd-d43dedba97cb
md"""
모든 임의의 ${\bf y}$에 대하여 아래가 성리하므로.

$${\bf y}^\top{\bf A}{\bf y}={\bf y}^\top{\bf \Psi}{\bf \Lambda}{\bf \Psi}^\top{\bf y}={\bf x}^\top {\bf \Lambda}{\bf x}=\sum_{i=1}^{n}x_i^2\lambda_i>0$$

단 여기에서 ${\bf x}={\bf \Psi}^\top{\bf y}$ 이다. 
"""

# ╔═╡ 350200af-c763-4bc3-8b7f-be8422988009
md"""
> 이전증명을 삭제하고 이 증명으로 수정하였습니다. 이 증명으로 공부하세요!!
"""

# ╔═╡ 706c2e23-626e-44e0-a443-70023ffbccc6
md"""
`-` 관찰2: 어떠한 실대칭행렬 ${\bf A}$가 positive definite matrix 이면 (즉 모든 non-zero vector ${\bf y}$에 대하여 ${\bf y}^\top {\bf A}{\bf y}>0$이 성립한다면) ${\bf A}$의 모든 고유값은 양수이다. 
"""

# ╔═╡ f8ccc269-0ea0-4601-ab4b-04b6aa1526ec
md"""
(왜?) 

행렬 ${\bf A}$는 실대칭행렬이므로 직교대각화가 가능하다. 일단 서로 직교하는 $n$개의 고유벡터 $\psi_1,\dots \psi_n$을 확보할 수 있다. $\psi_1$에 대하여 

$\psi_1^\top {\bf A}\psi_1=\psi_1^\top \lambda_1\psi_1=\lambda_1>0$

가 성립한다. 여기에서 첫번째 등호는 고유벡터의 정의, 두번째 등호는 $\Psi$가 직교행렬이라는 사실, 마지막의 부등호는 positive definite의 정의에 의하여 성립한다. 따라서 $\lambda_1>0$이다. 이러한 논의가 $\psi_2,\dots \psi_n$에 대하여 성립하므로 모든 고유값은 양수가 된다. 
"""

# ╔═╡ b44292d9-2a3f-49d6-bcc2-8fb5b72e9e6c
md"""
`-` 요약: ${\bf A}$가 실수인 대칭행렬이면 아래가 동치이다. 
- 행렬 ${\bf A}$의 모든 고유값이 양수이다. 
- 행렬 ${\bf A}$가 positive definite matrix이다. 
"""

# ╔═╡ da1a6561-ef07-492c-8e20-bfb0dc530a05
md"""
`-` 요약: ${\bf A}$가 실수인 대칭행렬이면 아래가 동치이다. 
- 행렬 ${\bf A}$의 모든 고유값이 0또는 양수이다. 
- 행렬 ${\bf A}$가 positive semidefinite matrix이다.
"""

# ╔═╡ a23df07e-24df-4ae3-942d-7ccfe024de57
md"""
## 숙제  

고유벡터행렬이 직교행렬이고 모든 고유값이 실수가 되는 매트릭스 ${\bf A}_{4\times 4}$ 를 구해보라. 

"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
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

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "0f4e115f6f34bbe43c19751c90a38b2f380637b9"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.3"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

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

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "1285416549ccfcdf0c50d4997a94331e88d68413"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╟─27d01165-2aa4-4b86-93ba-c7109636e692
# ╟─1d6ec17b-eca2-4289-8cfe-04a67bcac2c6
# ╟─d3ccc0c6-e092-11ec-3296-b5a32d51eade
# ╠═b868981c-dffb-4165-8e69-97889bc25312
# ╠═26300fa0-9f73-4874-9e60-d7b82f05cc71
# ╟─12aa1c0d-8f14-48cf-b244-c3fac884d676
# ╟─43a2c50f-0427-4c2c-af3e-c83d92c8e82c
# ╟─742c33c4-3bd5-4388-a9dd-4f9341eeb3e4
# ╠═9b52b105-40b4-4b77-b274-2d0a6914985a
# ╟─f1c60bf5-2052-4bf3-8e3b-ddb0cbd14b39
# ╠═0b882110-64b9-49c7-8ffd-66b95f0b8a5d
# ╟─bd476cf8-9946-4911-bee9-5170aabd9ca7
# ╠═2430ceb6-1154-4345-8487-05bd9a8136e2
# ╟─05af5915-30a6-49fd-850a-352a41036984
# ╠═9022b592-6ab0-40de-b96b-aecb56c663a5
# ╟─f080d7b0-4c1f-411e-9e86-d422dd77d2e7
# ╟─c226e8e1-e97c-4039-8baa-f82a00254917
# ╟─031acc24-0027-4db0-9270-b63a026e9575
# ╟─eac996bd-5458-40eb-a82d-88f3aac35033
# ╠═4d737613-177c-45e8-abe7-597ffd863049
# ╠═df004a39-17f0-47eb-9725-d46dd844f205
# ╠═4eabb492-2f64-4eb1-9041-a8d29eb40ff8
# ╟─5b94d293-7533-4f4d-8c71-0a2c19e94117
# ╟─1bcbd243-7088-4133-a8a9-6dab8d37001c
# ╟─94b5adf6-cbe5-4052-a2bd-dfabac6d7406
# ╠═6cc96118-41ac-4807-8f88-4ae86e491c79
# ╟─dc456627-72c5-4acf-89af-aa1ba01d2ba4
# ╟─fb480817-66f5-468b-b6c6-9bf145fc86c4
# ╟─3e069b39-28fb-4404-b292-8a54c2a25a8a
# ╟─7befd7ef-7ab9-4a00-a9d8-52cde80c23e6
# ╟─71d85e9f-a649-4584-8cbd-d43dedba97cb
# ╟─350200af-c763-4bc3-8b7f-be8422988009
# ╟─706c2e23-626e-44e0-a443-70023ffbccc6
# ╟─f8ccc269-0ea0-4601-ab4b-04b6aa1526ec
# ╟─b44292d9-2a3f-49d6-bcc2-8fb5b72e9e6c
# ╟─da1a6561-ef07-492c-8e20-bfb0dc530a05
# ╟─a23df07e-24df-4ae3-942d-7ccfe024de57
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
