### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 79ee3d48-46aa-42d2-b3b8-bf12bde667b6
using PlutoUI

# ╔═╡ e402f73a-264d-42cf-be51-b5fc45e05a1a
md"""
# 11wk-1: 벡터미분
"""

# ╔═╡ 6e5d5eda-1194-420c-81e2-d6ce8756331e
md"""
## 1. 강의영상
"""

# ╔═╡ 06b7582e-07de-4569-b2f7-75dfd50b365a
md"""
## 2. Imports
"""

# ╔═╡ 0da63a00-6cf6-47a4-a4c6-191a4df26837
PlutoUI.TableOfContents()

# ╔═╡ c5331296-1189-11ef-38cb-9bda1b37621c
md"""
## 3. 벡터미분
"""

# ╔═╡ 27a686c7-81d3-46a0-9309-5ec35254f55e
md"""
### A. 정의
"""

# ╔═╡ 69669497-8df3-42a3-b63a-c489aa077374
md"""
!!! info "벡터미분" 
	임의의 벡터 ${\bf y}_{n \times 1}$를 고려하자. 벡터미분은 아래와 같은 **벡터비슷한것**으로 정의한다. 
	
	$$\frac{\partial}{\partial \bf y}=\begin{bmatrix} \frac{\partial}{\partial y_1} \\ \frac{\partial}{\partial y_{2}} \\ \dots  \\ \frac{\partial}{\partial y_n} \end{bmatrix}$$

	대부분의 경우에서 $\frac{\partial}{\partial \bf y}$ 은 (n,1) col-vector 처럼 생각해도 무방하다. (그렇지만 col-vector 자체는 아니다) 
"""

# ╔═╡ 8f9b0d17-aa48-449d-a1ec-b0ee7c633550
md"""
!!! info "스칼라를 벡터로 미분" 
	임의의 스칼라 $z$와 벡터 ${\bf y}_{n\times 1}$을 고려하자. 스칼라를 벡터로 미분하는 경우 아래와 같은 방식으로 표현할 수 있다. 
	
	$$\frac{\partial}{\partial \bf y}  z := \begin{bmatrix} \frac{\partial z}{\partial y_1} \\ \frac{\partial z}{\partial y_{2}} \\ \dots  \\ \frac{\partial z}{\partial y_n} \end{bmatrix}$$

	스칼라에 대한 벡터미분은 언제나 잘 정의된다. 
"""

# ╔═╡ de35247b-88da-40a6-8e48-d997ba14f596
md"""
!!! info "벡터를 벡터로 미분" 
	
	임의의 벡터 ${\bf y}$와 ${\bf z}$를 고려하자. 벡터를 벡터로 미분하는 경우 아래와 같은 방식으로 표현 할 수 있다. 
	
	$$\frac{\partial }{\partial \bf y}\bf z^\top:=\frac{\partial \bf z^\top}{\partial \bf y}
	=\begin{bmatrix} 
	\frac{\partial z_1}{\partial y_1} & \frac{\partial z_2}{\partial y_1} & \dots  & \frac{\partial z_n}{\partial y_1} \\
	\frac{\partial z_1}{\partial y_2} & \frac{\partial z_2}{\partial y_2} & \dots  & \frac{\partial z_n}{\partial y_2} \\
	\dots & \dots & \dots  & \dots \\
	\frac{\partial z_1}{\partial y_n} & \frac{\partial z_2}{\partial y_n} & \dots  & \frac{\partial z_n}{\partial y_n}
	\end{bmatrix}$$

	벡터를 벡터로 미분하는 경우는 항상 잘 정의되는것은 아니다. 위의 미분 $\frac{\partial }{\partial \bf y}\bf z^\top$가 잘정의되려면 $\frac{\partial }{\partial \bf y}$를 벡터로 생각했을때 ${\bf z}^\top$와 행렬곱이 가능해야 한다. 위의 예시의 경우 $\frac{\partial }{\partial \bf y}$ 는 (n,1) 벡터로 생각할 수 있고, $\bf z^\top$ 은 (1,n) 벡터로 생각할 수 있어서 이들의 곱은 (n,n) 매트릭스로 생각할 수 있는데, 이때야 비로소 $\frac{\partial \bf z^\top}{\partial \bf y}$ 라는 표현을 쓸 수 있다. 
"""

# ╔═╡ 50e3404f-3368-45d1-a829-7169d6dfc92b
md"""
!!! warning "공식 억지로 외울필요 없음"
	언급하였듯이 그냥 $\frac{\partial}{\partial \bf y}$자체를 $n\times 1$ 매트릭스로 생각하고 $\bf z^\top$를 $1\times n$ 매트릭스로 생각하는 것이 속편하다. 즉 아래와 같이 생각하자. 
	
	$$\frac{\partial }{\partial \bf y}{\bf z}^\top=
	\begin{bmatrix}
	\frac{\partial}{\partial y_1}\\
	\dots\\
	\frac{\partial}{\partial y_n}\\
	\end{bmatrix}
	[z_1,\dots,z_n]=\begin{bmatrix} 
	\frac{\partial z_1}{\partial y_1} & \frac{\partial z_2}{\partial y_1} & \dots  & \frac{\partial z_n}{\partial y_1} \\
	\frac{\partial z_1}{\partial y_2} & \frac{\partial z_2}{\partial y_2} & \dots  & \frac{\partial z_n}{\partial y_2} \\
	\dots & \dots & \dots  & \dots \\
	\frac{\partial z_1}{\partial y_n} & \frac{\partial z_2}{\partial y_n} & \dots  & \frac{\partial z_n}{\partial y_n}
	\end{bmatrix}$$

	이러면 공식 안외워도 된다.
"""

# ╔═╡ f80913c3-86bf-4b71-b99d-d25867fae118
md"""
### B. 연습문제
"""

# ╔═╡ be4e1e85-e091-4ae7-974d-772f8a966d17
md"""
-- 예제1: 아래의 미분들의 계산결과를 보여라. 여기에서 각 벡터 및 매트릭스의 차원은 ${\bf x}_{n\times 1}, {\bf y}_{n\times 1}, {\bf X}_{n\times p}, {\boldsymbol \beta}_{p\times 1}$ 로 가정한다. 

`(1)` ``\frac{\partial}{\partial {\bf x}}{\bf x}^\top {\bf y} = {\bf y}``

`(2)` ``\frac{\partial}{\partial {\bf x}}{\bf y}^\top {\bf x} = {\bf y}``

`(3)` ``\frac{\partial}{\partial {\boldsymbol \beta}}{\bf y}^\top {\bf X} {\boldsymbol \beta} = {\bf X}^\top{\bf y}``
"""

# ╔═╡ 247fd51b-7049-4da5-8fde-599754da268f
md"""
-- 예제2: 아래 미분들의 계산결과를 보여라 .여기에서 각 벡터 및 매트릭스의 차원은 ${\bf y}_{n\times 1}, {\bf X}_{n\times p}, {\boldsymbol \beta}_{p\times 1}$ 로 가정한다. 

`(1)` $\frac{\partial}{\partial {\bf y}}{\bf y}^\top{\bf y}=2{\bf y}$

`(2)` $\frac{\partial}{\partial {\boldsymbol \beta}}{\boldsymbol \beta}^\top{\bf X}^\top{\bf X}{\boldsymbol \beta}=2{\bf X}^\top{\bf X}{\boldsymbol \beta}$
"""

# ╔═╡ 2e788885-a662-4567-9406-7e2e5f9de351
md"""
-- 예제3: 아래와 같은 함수를 최소화하는 ${\boldsymbol \beta}$를 구하라. 여기에서 각 벡터 및 매트릭스의 차원은 ${\bf y}_{n\times 1}, {\bf X}_{n\times p}, {\boldsymbol \beta}_{p\times 1}$ 로 가정한다. 

$$loss:=({\bf y}-{\bf X}{\boldsymbol \beta})^\top({\bf y}-{\bf X}{\boldsymbol \beta})$$

"""

# ╔═╡ f65ed75c-fbc1-4b47-a5de-ca0ac1a67781
md"""
-- 예제4: 아래와 같은 함수를 최소화하는 ${\boldsymbol \beta}$를 구하라. 여기에서 각 벡터 및 매트릭스의 차원은 ${\bf y}_{n\times 1}, {\bf X}_{n\times p}, {\boldsymbol \beta}_{p\times 1}$ 로 가정한다. 

$$loss:=({\bf y}-{\bf X}{\boldsymbol \beta})^\top({\bf y}-{\bf X}{\boldsymbol \beta})+\lambda{\boldsymbol \beta}^\top{\boldsymbol \beta}$$
"""

# ╔═╡ e7048c45-905b-43de-a11a-2f21dddd61a2
md"""
-- 예제5: 예제3과 같은 함수를 아래와 같이 재표현하자. 

- ``{\boldsymbol u} = {\bf X}{\boldsymbol \beta}``
- ``{\boldsymbol v} = {\bf y}-{\boldsymbol u}``
- ``loss = {\boldsymbol v}^\top {\boldsymbol v}``

아래와 같이 주장할 수 있는가?

$$\frac{\partial}{\partial \boldsymbol \beta} loss =\bigg(\frac{\partial}{\partial \boldsymbol \beta}{\boldsymbol u}^\top\bigg)\bigg(\frac{\partial}{\partial \boldsymbol u}{\boldsymbol v}^\top\bigg)\bigg(\frac{\partial}{\partial \boldsymbol v}loss\bigg)$$

"""

# ╔═╡ 4bfb8c49-9978-4c6e-b2b7-f2576e9e98ed
md"""
-- 예제6: (직교행렬이 아닌) 임의의 ${\bf V}$ 와 임의의 대각행렬 ${\bf \Lambda}$에 대하여 아래를 보여라. 여기에서 각 매트릭스의 차원은 ${\bf X}_{n\times p}, {\bf V}_{p\times p}, {\bf \Lambda}_{p \times p}$ 로 가정한다. 

`(1)` ``\frac{\partial}{\partial {\bf V}}tr({\bf V}^\top{\bf X}^\top{\bf X}{\bf V}) = 2{\bf X}^\top{\bf X}{\bf V}``

`(2)` ``\frac{\partial}{\partial {\bf V}}tr({\bf \Lambda}({\bf V}^\top{\bf V}-{\bf I})) = 2{\bf V}{\bf \Lambda}``

`(3)` ``\frac{\partial}{\partial {\bf \Lambda}}tr({\bf \Lambda}({\bf V}^\top {\bf V}-{\bf I}))={\bf V}^\top{\bf V}-{\bf I}``
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.3"
manifest_format = "2.0"
project_hash = "6e7bcec4be6e95d1f85627422d78f10c0391f199"

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

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

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

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

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
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

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

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─e402f73a-264d-42cf-be51-b5fc45e05a1a
# ╟─6e5d5eda-1194-420c-81e2-d6ce8756331e
# ╟─06b7582e-07de-4569-b2f7-75dfd50b365a
# ╠═79ee3d48-46aa-42d2-b3b8-bf12bde667b6
# ╠═0da63a00-6cf6-47a4-a4c6-191a4df26837
# ╟─c5331296-1189-11ef-38cb-9bda1b37621c
# ╟─27a686c7-81d3-46a0-9309-5ec35254f55e
# ╟─69669497-8df3-42a3-b63a-c489aa077374
# ╟─8f9b0d17-aa48-449d-a1ec-b0ee7c633550
# ╟─de35247b-88da-40a6-8e48-d997ba14f596
# ╟─50e3404f-3368-45d1-a829-7169d6dfc92b
# ╟─f80913c3-86bf-4b71-b99d-d25867fae118
# ╟─be4e1e85-e091-4ae7-974d-772f8a966d17
# ╟─247fd51b-7049-4da5-8fde-599754da268f
# ╟─2e788885-a662-4567-9406-7e2e5f9de351
# ╟─f65ed75c-fbc1-4b47-a5de-ca0ac1a67781
# ╟─e7048c45-905b-43de-a11a-2f21dddd61a2
# ╟─4bfb8c49-9978-4c6e-b2b7-f2576e9e98ed
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
