### A Pluto.jl notebook ###
# v0.19.4

using Markdown
using InteractiveUtils

# ╔═╡ 936acea6-fc7f-4e12-b9fe-5b354fe08f9d
using LinearAlgebra, PlutoUI

# ╔═╡ 11333c04-915c-49df-8ac5-d85bed1d9050
md"""
# 6월2일
"""

# ╔═╡ 4c601d48-0c6b-4405-b002-e3868417373a
html"""
<div style="display: flex; justify-content: center;">
<div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
<iframe src=
"
https://www.youtube.com/embed/playlist?list=PLQqh36zP38-wyEPkRBxbkVcqY5I4HRU1u
"
width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ 471d48be-e20a-11ec-37b6-51839c9abe25
md"""
## usings
"""

# ╔═╡ 83a09e08-851e-49a4-a395-e81ef92fdedd
PlutoUI.TableOfContents()

# ╔═╡ 05d65e58-f80e-41e1-b38b-8bee321947e4
md"""
## SVD
"""

# ╔═╡ 1b89725b-a3fc-4961-9f36-2d63b1cf22ba
md"""
### SDV의 손계산
`-` ${\bf X}={\bf U}{\bf D}{\bf V}^\top$를 계산하는 방법 (손으로 하는 방법)
- step1: ${\bf X}^\top {\bf X}$의 고유벡터행렬을 구하고 그것을 이용하여 ${\bf V}$를 구한다. 
- step2: ${\bf X}{\bf X}^\top$의 고유벡터행렬을 구하고 그것을 이용하여 ${\bf U}$를 구한다.
- step3: ${\bf X}^\top{\bf X}$의 0이 아닌 고유값 혹은 ${\bf X}{\bf X}^\top$의 0이 아닌 고유값을 구하고 (두 행렬의 0이 아닌 고유값은 일치함) 루트를 취한다. 고유값을 큰 순서대로 나열하여 ${\bf D}$ matrix의 대각선에 넣는다. 
"""

# ╔═╡ e8fa1778-1661-467e-bd7b-53c63fcacdbc
md"""
`-` 왜 위와 같은 step으로 계산하면 svd가 되는지는 선형대수학 시간에.. 
"""

# ╔═╡ e1fbbe71-0efa-4e1a-9425-fbe78724b1cd
md"""
`-` 아래의 포인트를 파악하는 것이 중요. 
- 포인트1: ${\bf X}$의 v-matrix (right singular matrix) 는 ${\bf X}^\top {\bf X}$의 고유벡터 행렬에서 구할 수 있고 ${\bf X}$의 u-matrix (left singular matrix) 는 ${\bf X}{\bf X}^\top$의 고유벡터 행렬에서 구할 수 있다. // (이때 ${\bf X}^\top {\bf X}$ 혹은 ${\bf X}{\bf X}^\top$ 는 실대칭행렬이므로 모든 column들이 서로 직교하고 full rank인 고유벡터행렬을 가질 수 있다, recall the spectral thm.)
- 포인트2: ${\bf X}$의 d-matrix는 모든 대각선의 원소가 non-negative 이다. 왜냐하면 ${\bf X}^\top{\bf X}$와 ${\bf XX}^\top$의 0이 아닌 고유값을 구하고 루트를 취하여 크기가 큰 순서대로 나열했기 때문에. (이때 당당하게 루트를 취할 수 있는 이유는 두 매트릭스 모두 real, symm, psd matrix 이기 때문!)
"""

# ╔═╡ 051f2a9e-fab3-400b-8ef8-b6fe0ae73e02
md"""
### real-symm and psd matrix의 SVD 
"""

# ╔═╡ ea557115-904f-4d28-a136-8f700d3ffb41
md"""
`-` 이론: (정사각)행렬 ${\bf A}_{n\times n}$가 실대칭행렬이고 positive semidefinite matrix 이면 ${\bf A}$의 고유값분해와 SVD가 일치한다. 즉 아래를 만족하는 ${\bf \Psi}, {\bf \Lambda}, {\bf U}, {\bf D}, {\bf V}$ 를 선택할 수 있다. 

${\bf A}={\bf \Psi}{\bf \Lambda}{\bf \Psi}^\top={\bf U}{\bf D}{\bf V}^\top$

-  ${\bf \Psi}={\bf U}={\bf V}$ 이고 이들은 모두 직교행렬. 
-  ${\bf \Lambda}={\bf D}$ 이고 이들은 모두 대각선의 원소가 non-negative인 대각행렬.
- note: 모든 매트릭스는 $n\times n$ 매트릭스. 
"""

# ╔═╡ df76cae4-d575-440c-8d78-22b68af23a1d
md"""
`-` 아래를 따지면서 위의 이론을 체크해보자. 

(1) ${\bf A}$의 u-matrix와 v-matrix는 같다. 

(2) ${\bf A}$의 고유벡터행렬은 ${\bf A}$ v-matrix의 정의를 만족한다. 

(3) ${\bf A}$의 v-matrix는 ${\bf A}$의 고유벡터행렬의 정의를 만족한다. 
"""

# ╔═╡ 1c0929a1-4cb7-4a3f-91dd-39214d28579f
md"""
### 예제풀이
"""

# ╔═╡ 84941afa-70cc-4884-958c-d0926535de6e
X = [1 2
	 3 4
	 5 6]

# ╔═╡ 65cfa705-b51c-48a9-b2b2-1c91e64fff3e
X'X # 이렇게 계산해도 곱셈이 수행된다. 

# ╔═╡ 412f387f-5632-4752-90d6-667ac0edc9ae
md"""
(예제1) ${\bf X}$의 svd를 이용하여 ${\bf X}^\top{\bf X}$의 고유값분해를 구하라.
"""

# ╔═╡ d242395d-d819-4793-bcb3-1fc07662b965
let 
	U,d,V = svd(X)

	Ψ = V 
	Λ = Diagonal(d.^2)
	Ψ * Λ * Ψ'
end 

# ╔═╡ 1942a3bb-0403-4aac-b15e-5d30ed677c52
md"""
(예제2) ${\bf X}^\top {\bf X}$의 고유값분해를 이용하여 ${\bf X}^\top {\bf X}$의 svd를 구하라. 
"""

# ╔═╡ 79f5a917-6c17-4c2d-a28a-cd06499d2dd5
let 
	λ, Ψ = eigen(X'X)

	U = Ψ
	V = Ψ 
	D = Diagonal(λ)
	U*D*V'
end

# ╔═╡ 5e8238b2-91a6-4f21-a503-c8d7ccfa08d4
md"""
(예제3) ${\bf X}^\top{\bf X}$의 고유값분해와 ${\bf X}{\bf X}^\top$의 고유값분해를 이용하여 ${\bf X}$의 svd를 구하라.
"""

# ╔═╡ 1a306a2e-6bf1-43da-bef8-ca426c18a139
let 
	λu, Ψu = eigen(X*X') # Ψu is 3*3 matrix , λu is length 3 vector 
	λv, Ψv = eigen(X'X) # Ψv is 2*2 matrix , λv is length 2 vector 
	## svd의 ver1로 U,D,V matrix를 구성하면 편리할 것 같다. 
	d= [√(λv[2]), √(λv[1])] # 0이 아닌 고유값을 큰 순서대로 나열하고 루트를 취함
	D = [Diagonal(d); [0 0]] # D는 3*2 matrix이어야함
	U = [Ψu[:,3] Ψu[:,2] Ψu[:,1]] # 고유값을 큰 순서대로 나열했으니 벡터도 재정렬
	V = [Ψv[:,2] -Ψv[:,1]] # 고유값을 큰 순서대로 나열했으니 고유벡터도 재정렬
	U * D * V' 
end 

# ╔═╡ f3fd2f61-c458-4954-ac3d-714d7335775c
md"""
- 교훈: 모든 예제가 이론만큼 코딩이 깨끗하지는 않다. 
- 그래서 예제1,2 정도만 잘 기억하면 된다. 
"""

# ╔═╡ 63283701-f861-4ff3-8982-05a9019f29d6
md"""
## 숙제
"""

# ╔═╡ 85c39a5d-df42-4c28-a7da-cd41fab5dde9
md"""
다음을 잘 읽고 참 거짓을 판단하여라. (단, 매트릭스의 모든 원소는 실수라고 가정한다.)
"""

# ╔═╡ b8b90462-4bc7-4b42-9e1b-9031a32f3fa8
md"""
(1) ${\bf A}$가 실대칭행렬이면 ${\bf A}$의 모든 고유값은 실수이다. 

(2) ${\bf A}$가 실대칭행렬이면 ${\bf A}$의 고유벡터행렬을 직교행렬로 선택할 수 있다. 

(3) 실대칭행렬 ${\bf A}$의 모든 고유값이 양수이면 ${\bf A}$는 positive definite matrix 이다. 

(4) 실대칭행렬 ${\bf A}$의 고유벡터행렬은 항상 full rank 이다. 

(5) 임의의 행렬 ${\bf X}$에 대하여 ${\bf X}^\top{\bf X}$는 항상 대각화가능하다. 
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
# ╟─11333c04-915c-49df-8ac5-d85bed1d9050
# ╟─4c601d48-0c6b-4405-b002-e3868417373a
# ╟─471d48be-e20a-11ec-37b6-51839c9abe25
# ╠═936acea6-fc7f-4e12-b9fe-5b354fe08f9d
# ╠═83a09e08-851e-49a4-a395-e81ef92fdedd
# ╟─05d65e58-f80e-41e1-b38b-8bee321947e4
# ╟─1b89725b-a3fc-4961-9f36-2d63b1cf22ba
# ╟─e8fa1778-1661-467e-bd7b-53c63fcacdbc
# ╟─e1fbbe71-0efa-4e1a-9425-fbe78724b1cd
# ╟─051f2a9e-fab3-400b-8ef8-b6fe0ae73e02
# ╟─ea557115-904f-4d28-a136-8f700d3ffb41
# ╟─df76cae4-d575-440c-8d78-22b68af23a1d
# ╟─1c0929a1-4cb7-4a3f-91dd-39214d28579f
# ╠═84941afa-70cc-4884-958c-d0926535de6e
# ╠═65cfa705-b51c-48a9-b2b2-1c91e64fff3e
# ╟─412f387f-5632-4752-90d6-667ac0edc9ae
# ╠═d242395d-d819-4793-bcb3-1fc07662b965
# ╟─1942a3bb-0403-4aac-b15e-5d30ed677c52
# ╠═79f5a917-6c17-4c2d-a28a-cd06499d2dd5
# ╟─5e8238b2-91a6-4f21-a503-c8d7ccfa08d4
# ╠═1a306a2e-6bf1-43da-bef8-ca426c18a139
# ╟─f3fd2f61-c458-4954-ac3d-714d7335775c
# ╟─63283701-f861-4ff3-8982-05a9019f29d6
# ╟─85c39a5d-df42-4c28-a7da-cd41fab5dde9
# ╟─b8b90462-4bc7-4b42-9e1b-9031a32f3fa8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
