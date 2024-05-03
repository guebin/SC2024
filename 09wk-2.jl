### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 2db4df04-0877-42f7-a2f4-1e7d82da88ea
using LinearAlgebra, PlutoUI

# ╔═╡ 73aa567a-082c-11ef-19ff-8947f8f86252
md"""
# 09wk-2: 실대칭행렬, PD-행렬, SVD의 존재
"""

# ╔═╡ dc9a32c9-366b-453b-a876-d13f8b7d29f3
md"""
## 1. 강의영상
"""

# ╔═╡ eb84030e-e9c1-48c7-8f9e-138b7e4adfee
# html"""
# <div style="display: flex; justify-content: center;">
# <div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
# <iframe src=
# "
# https://youtube.com/embed/playlist?list=PLQqh36zP38-y-PZZgrJ4qHRKk82p82Hj4&si=YxiCdULCVSmT_lfD
# "
# width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
# """

# ╔═╡ f9ccb632-23d6-4531-8e3f-b076685baeb6
md"""
## 2. Imports
"""

# ╔═╡ 586aeb22-c4d3-4eb0-8ed7-49ccebb625a1
PlutoUI.TableOfContents()

# ╔═╡ f32588dc-ebeb-4477-8225-045154104adb
md"""
## 3. 실대칭행렬
"""

# ╔═╡ e7c03452-fb03-4907-b689-8e2afd587961
md"""
### A. 스펙트럼 정리
"""

# ╔═╡ d0731163-6d63-4abd-a1db-11e6b462559a
md"""
!!! info "이론: 스펙트럼 정리 (Spectral Theorem)" 
	임의의 정사각행렬 ${\bf A}$에 대하여 아래가 성립한다. 

	- ``{\bf A}`` 가 실대칭행렬 $\Leftrightarrow$ ${\bf A}$는 (1) 모든 고유값이 실수이고 (2) 직교대각화 가능함. 

	여기에서 $\Rightarrow$ 방향을 스펙트럼정리라고 한다. (반대방향은 당연해서..) 스펙트럼 정리를 엄밀하게 다시쓰면 아래와 같다. 

	- ``{\bf A}`` 가 실대칭행렬 $\Leftrightarrow$  ${\bf A}={\bf \Psi} {\bf \Lambda}{\bf \Psi}^\top$ 를 만족하는 실대각행렬 ${\bf }$를 $$\bf \Psi$와 

	즉 임의의 ``n \times n`` 행렬 ``{\bf A}``에 대하여 아래가  실대칭행렬이라면 ${\bf A}$를 아래와 같이 표현할 수 있으며 역도 성립한다. 

	${\bf A} = {\bf \Psi} {\bf \Lambda}{\bf \Psi}^\top$

	여기에서 ${\bf \Psi}$는 ${\bf A}$의 고유벡터행렬이고, ${\bf \Lambda}$는 고유값행렬이다. 
"""

# ╔═╡ 0d3ae451-7512-4016-8394-843091a71767
md"""
### B. 응용
"""

# ╔═╡ 62c5757a-6807-4c73-9b7f-88128b0d4370
md"""
-- 예제1. 실대칭행렬일때 "고유값분해 = 특이값분해" 가 성립하지 않는 이유를 설명해보라. 
"""

# ╔═╡ 8df6b3c7-6563-475d-88e7-19100f30488e
md"""
-- 예제2. 아래의 행렬들이 (1) 모든 고유값이 실수인지 (2) 직교대각화 가능한지 체크하라.
"""

# ╔═╡ b8ec122e-c5a8-42c2-b579-d7e12ff62a5a
let 
	A = [1 3
		 3 2]
end 

# ╔═╡ 4ee89af2-e782-4e25-aba8-40841d5431aa
let 
	A = [1 -im
		 im  1]
	λ,Ψ = eigen(A)
end 

# ╔═╡ d6aed637-0f86-46ee-94bb-0b0225c2bb74


# ╔═╡ 7b6c2413-79c7-491c-b872-d15fa9944685
md"""
-- 예제2. 직교대각화 가능한 행렬 ${\bf A}$의 고유값중 0이 하나도 없다면, ${\bf A}$의 역행렬이 존재함을 보여라. 
"""

# ╔═╡ 470c3022-480f-4291-9b9b-0ff8f6efe694
# 항상 아래와 같이 찾을 수 있음 

# ╔═╡ ee73c7f3-20d7-406b-8279-80b4268443d4
let 
	A = 

end 

# ╔═╡ 1b7fe945-7b48-42d0-92cc-1130d7fb0d1f
md"""
-- 예제3. 예제2의 논리가 대각화가능행렬 ${\bf A}$에 대하여서도 동일하게 성립함을 보여라. 즉 대각화가능행렬 ${\bf A}$의 모든 고유값이 0이 아니라면, ${\bf A}$의 역행렬이 존재함을 보여라. 
"""

# ╔═╡ d480d5f3-53ae-4fb6-877c-60206626f80f
md"""
## 4. PD행렬
"""

# ╔═╡ 17fc6337-1ad4-40bf-8887-9930ec76a21e
md"""
### A. PD/SPD, ND/NSD matrix의 정의
"""

# ╔═╡ 06f5d756-8656-431d-834d-b6f591fe4967
md"""
!!! info "정의 : PD/PSD, ND/NSD" 
	임의의 non-zero real-vector ${\bf y}_{n\times 1}$ 에 대하여 아래를 만족하는 실대칭행렬 ${\bf A}_{n \times n}$를 positive definite matrix 라고 한다. 

	${\bf y}^\top {\bf A} {\bf y} > 0$ 

	이외에도 $\geq, < ,\leq$ 에 따라 postive semi-definite, negative definite, negative semi-definite matrix를 정의할 수 있다. 
"""

# ╔═╡ 1d247622-5ab9-4207-9b49-35e5c0eae9e4
md"""
### B. PD-행렬의 성질
"""

# ╔═╡ 21588176-bc55-4d79-82b2-6192e3871ffc
md"""
!!! info "정리: " 
	행렬 ``{\bf A}_{n\times n}`` 가 PD-행렬이라면 (1) 모든 고유값이 양수이고 (2) 직교대각화 가능하다. 그리고 ${\bf A}_{n\times n}$ 이 PSD, ND, NSD 에 대하여서도 비슷한 정리가 성립한다. 
"""

# ╔═╡ f6e131e4-d5a6-488c-bb71-78aa9ab09425
md"""
### C. 생각해 볼 점
"""

# ╔═╡ f25cb325-a290-4156-b978-c4de13840b20
md"""
-- 예제1: 임의의 행렬 ${\bf X}_{n\times p}$에 대하여 ${\bf X}^\top{\bf X}$ 와 ${\bf X}{\bf X}^\top$는 PD-행렬임을 보여라.
"""

# ╔═╡ 42a677cc-06b4-473c-be9d-00ecee71448e
md"""
-- 예제2: 임의의 행렬 ${\bf A}$가 PD-행렬이라면 아래를 만족하는 행렬들이 존재함을 보여라. (찾으면 됩니당)

- ``{\bf A}^{1/2}{\bf A}^{1/2}={\bf A}``
- ``{\bf A}^{-1/2}{\bf A}^{-1/2}={\bf A}^{-1}``
"""

# ╔═╡ e7fd9ae0-d22e-4a77-8a9e-3e1969f18022
md"""
-- 예제3: 아래가 PD-행렬임을 보이고 ${\bf A}^{-1/2}$와 ${\bf A}^{1/2}$을 구하라. 
"""

# ╔═╡ c5525c58-90c4-4a1a-b038-2abd291a5e10
let 
	A = [1.0 0.5 
		 0.5 1.0]
end

# ╔═╡ 1186d3d4-c5cc-4eeb-8d44-eae2aaabd43a
md"""
## 5. Story 정리 -- 대충정리
"""

# ╔═╡ 3d4fba9d-5e85-4f1d-8fc8-865b5f839d2f
md"""
`1`. ${\bf A}$가 정사각행렬이면, ${\bf A}{\bf \Psi}={\bf \Psi}{\bf \Lambda}$ 로 쓸 수 있음. 

`2`. ${\bf A}$가 정사각행렬이고 ${\bf A}$의 고유벡터행렬이 full-rank 이면, ${\bf A}={\bf \Psi}{\bf \Lambda}{\bf \Psi}^{-1}$ 라고 쓸 수 있음. 

`3`. ${\bf A}$가 실대칭행렬이고 (따라서 정사각행렬이라는 소리)
"""

# ╔═╡ 73512dff-60dc-47b2-bac6-7d2faebb3656
md"""
## 6. SVD의 존재성 ($\star$)
"""

# ╔═╡ 3a60d7d9-e0be-40a6-8a8c-d6f766888a30
md"""
### A. SVD의 존재성
"""

# ╔═╡ 2d15f25e-bfe1-4787-830a-f0633e98828b
md"""
이말하려고 먼길 돌아왔습니다. 임의의 행렬 ${\bf X}$에 대하여 SVD는 무조건 존재합니다. 
"""

# ╔═╡ 26059dc2-ecb9-4f8e-845d-a4e0cee2bbc9
md"""
!!! info "정리: 특이값분해 (SVD)" 
	임의의 행렬 ``{\bf X}_{n\times p}`` 에 대하여 특이값 분해는 항상 존재한다. 
"""

# ╔═╡ 6c4d6e46-a19a-4e7c-9d42-361bcfad21bb
md"""
왜?
"""

# ╔═╡ 12ae53ef-924c-42d3-8a0a-28b221e2b726
md"""
### B. ${\bf U}$, ${\bf V}$의 정체
"""

# ╔═╡ 75776417-3359-4e3a-8f28-76d7b0a316a6
md"""
-- 예제1: ${\bf V}$는 ${\bf X}^\top{\bf X}$의 고유벡터임을 보여라.
"""

# ╔═╡ 68f8a231-03bf-47cb-a86e-3b4457a07259
md"""
-- 예제2: ${\bf U}$는 ${\bf X}{\bf X}^\top$의 고유벡터임을 보여라.
"""

# ╔═╡ 99428e3c-1e8f-408e-bcb1-0b266314becf
md"""
## 6. 특이값분해 = 고유값분해
"""

# ╔═╡ 7f04c611-6835-478d-a4a0-244f63837f2f
md"""
!!! info "정리: 특이값분해 = 고유값분해" 
	행렬 ``{\bf A}_{n\times n}`` 가 PSD-행렬이라면 특이값분해와 고유값분해가 일치하도록 만들어주는 
	- ``{\bf U}={\bf V}={\bf \Psi}``
	- ``{\bf D}={\bf \Lambda}``
	를 선택할 수 있다.
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.3"
manifest_format = "2.0"
project_hash = "e41793cbd1124ea5d05573eb874098b20a960e1d"

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
# ╟─73aa567a-082c-11ef-19ff-8947f8f86252
# ╟─dc9a32c9-366b-453b-a876-d13f8b7d29f3
# ╠═eb84030e-e9c1-48c7-8f9e-138b7e4adfee
# ╟─f9ccb632-23d6-4531-8e3f-b076685baeb6
# ╠═2db4df04-0877-42f7-a2f4-1e7d82da88ea
# ╠═586aeb22-c4d3-4eb0-8ed7-49ccebb625a1
# ╟─f32588dc-ebeb-4477-8225-045154104adb
# ╟─e7c03452-fb03-4907-b689-8e2afd587961
# ╠═d0731163-6d63-4abd-a1db-11e6b462559a
# ╟─0d3ae451-7512-4016-8394-843091a71767
# ╟─62c5757a-6807-4c73-9b7f-88128b0d4370
# ╟─8df6b3c7-6563-475d-88e7-19100f30488e
# ╠═b8ec122e-c5a8-42c2-b579-d7e12ff62a5a
# ╠═4ee89af2-e782-4e25-aba8-40841d5431aa
# ╠═d6aed637-0f86-46ee-94bb-0b0225c2bb74
# ╟─7b6c2413-79c7-491c-b872-d15fa9944685
# ╠═470c3022-480f-4291-9b9b-0ff8f6efe694
# ╠═ee73c7f3-20d7-406b-8279-80b4268443d4
# ╟─1b7fe945-7b48-42d0-92cc-1130d7fb0d1f
# ╟─d480d5f3-53ae-4fb6-877c-60206626f80f
# ╟─17fc6337-1ad4-40bf-8887-9930ec76a21e
# ╟─06f5d756-8656-431d-834d-b6f591fe4967
# ╟─1d247622-5ab9-4207-9b49-35e5c0eae9e4
# ╠═21588176-bc55-4d79-82b2-6192e3871ffc
# ╟─f6e131e4-d5a6-488c-bb71-78aa9ab09425
# ╟─f25cb325-a290-4156-b978-c4de13840b20
# ╟─42a677cc-06b4-473c-be9d-00ecee71448e
# ╟─e7fd9ae0-d22e-4a77-8a9e-3e1969f18022
# ╠═c5525c58-90c4-4a1a-b038-2abd291a5e10
# ╟─1186d3d4-c5cc-4eeb-8d44-eae2aaabd43a
# ╠═3d4fba9d-5e85-4f1d-8fc8-865b5f839d2f
# ╟─73512dff-60dc-47b2-bac6-7d2faebb3656
# ╟─3a60d7d9-e0be-40a6-8a8c-d6f766888a30
# ╟─2d15f25e-bfe1-4787-830a-f0633e98828b
# ╟─26059dc2-ecb9-4f8e-845d-a4e0cee2bbc9
# ╟─6c4d6e46-a19a-4e7c-9d42-361bcfad21bb
# ╟─12ae53ef-924c-42d3-8a0a-28b221e2b726
# ╟─75776417-3359-4e3a-8f28-76d7b0a316a6
# ╟─68f8a231-03bf-47cb-a86e-3b4457a07259
# ╟─99428e3c-1e8f-408e-bcb1-0b266314becf
# ╟─7f04c611-6835-478d-a4a0-244f63837f2f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
