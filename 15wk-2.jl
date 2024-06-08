### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 771b13aa-2bdd-400a-87e3-cc5bae383cfe
md"""
# 15wk-2: 기말고사
"""

# ╔═╡ 3a89c770-a071-4d46-bff3-57e54008f6e8
md"""
## 1. SVD
"""

# ╔═╡ 8f8f60db-9082-47a5-8de0-60f11572afcd
md"""
통계학에서 SVD는 어떻게 활용될 수 있는가? 활융분야를 목록화하고 간단히 서술하라. 
"""

# ╔═╡ aa306e3f-6e9a-413c-8e75-22f49fbc150a
md"""
## 2. 다중공선성
"""

# ╔═╡ 9bf50bb3-350c-4edd-9a9f-a41c0fc865c2
md"""
토익,텝스,학점을 설명변수 ${\bf X}_{n \times 3}$ 로 설정하고 이를 바탕으로 연봉 ${\bf y}_{n\times 1}$를 추정하고자 한다. (이때 학생들의 토익과 텝스점수는 서로 비슷하다고 가정한다. 즉 토익점수가 높은 학생은 대체로 텝스점수도 높으며, 반대의 경우도 그러하다고 가정한다) 다음을 잘 읽고 물음에 답하라. 
"""

# ╔═╡ 51e8f9e3-b7d4-4e78-8b92-a26af54cda8b
md"""
(1) 선형회귀를 사용하여 계수(토익,텝스,학점이 연봉에 미치는 영향)를 추정하고자 한다. 이러한 상황은 그림1에서 무엇과 관련이 있는가? 왜 그렇다고 생각하는가? 
"""

# ╔═╡ 6d14b11f-0f77-470a-9512-63fad98dd79a
md"""
![](https://www.kdnuggets.com/wp-content/uploads/arya_biasvariance_tradeoff_5.png)
ref: <https://www.kdnuggets.com/2022/08/biasvariance-tradeoff.html>
"""

# ╔═╡ e069c354-8ee3-4490-8937-25a665f0ce54
md"""
(2) 능형회귀를 이용하여 계수를 추정한다고 하자. 여기에서 $\lambda$는 어떠한 역할을 하는가? 그림과 연관시켜 설명하라.
"""

# ╔═╡ ad99f90c-512f-450e-89e5-3283b6d68a57
md"""
(3) 주성분 회귀 (Principal component regression, PCR)을 이용하여 계수를 추정하고자 한다. 이때 principle componet 수를 작게 설정할때와 크게 설정할때 어떠한 일이 생기는지 설명하라.
"""

# ╔═╡ 734e8d32-55de-4eaf-8498-10c9e9f2e51f
md"""
(4) 능형회귀에서 $\lambda=0$ 으로 설정하거나 $\lambda = \infty$로 설정하는 것이 어떠한 의미를 가지는 주성분 회귀와 연결시켜 설명하라.
"""

# ╔═╡ 14c22744-c549-4fac-84e0-7f913f1cfa58
md"""
## 3. 면접질문?
"""

# ╔═╡ ddd0c1bd-19fe-45ab-893b-c94242d612c0
md"""
(1) 능형회귀에 대하여 간단히 설명하라. 
"""

# ╔═╡ a329cfda-d89c-49e0-9f0e-59d2cad465d5
md"""
(2) 다중공선성이란 무엇이며 어떤 문제를 일으키는 간단히 서술하라.
"""

# ╔═╡ 88d2172a-9df3-454d-ae8b-290a069a6d27
md"""
(3) ${\bf X}_{n\times p}$, $p>2$ 일 경우 ${\bf X}$를 시각화하는 방법에 대하여 간단히 서술하라.
"""

# ╔═╡ 5d1a6abd-3075-4b7f-95a8-2a4101a81185
md"""
(4) 직교변환이 가지는 의미를 간단히 서술하라.
"""

# ╔═╡ e775138e-2889-4d3b-b9c8-bade656ebc43
md"""
(5) ``{\bf X}``가 이변량 정규분포를 따른다고 가정하자. $\mathbb{V}({\bf X})$의 고유벡터행렬을 활용하는 통계적 처리기법을 있는가? 있다면 서술하라. (하나만 서술해도 무방)
"""

# ╔═╡ 5232382a-6dd6-42e2-9e9b-92da7dbad672
md"""
(6) SVD를 이용하여 이미지를 압축하는 방법을 간단히 서술하라.
"""

# ╔═╡ 65e570f3-3597-449c-8e8e-7fe563cef9b9
md"""
(7) 주성분분석을 하게 되면 얻게되는 이점을 간단히 서술하라.
"""

# ╔═╡ ac002d2e-1064-433e-b695-01702083b6a4
md"""
(8) 선형변환을 SVD를 이용하여 해석하라.
"""

# ╔═╡ 1cd6bfe2-a720-46a3-adfe-bebfc23189ef
md"""
(9) 변환을 의미하는 행렬 ${\bf A}$가 데이터를 의미하는 행렬 ${\bf X}$의 앞에 곱해지는 경우와 뒤에 곱해지는 경우 각각 어떠한 의미를 가지는지 설명하라.
"""

# ╔═╡ 816f1e62-3b24-459a-86fb-bfa6624cbe57
md"""
(10) R(*lm()*)과 Python(*sklearn.linear_model*)에서 더미변수가 포함된 회귀분석을 수행하는 로직이 다르다. 차이점에 대하여 서술하라.
"""

# ╔═╡ Cell order:
# ╟─771b13aa-2bdd-400a-87e3-cc5bae383cfe
# ╟─3a89c770-a071-4d46-bff3-57e54008f6e8
# ╟─8f8f60db-9082-47a5-8de0-60f11572afcd
# ╟─aa306e3f-6e9a-413c-8e75-22f49fbc150a
# ╟─9bf50bb3-350c-4edd-9a9f-a41c0fc865c2
# ╟─51e8f9e3-b7d4-4e78-8b92-a26af54cda8b
# ╟─6d14b11f-0f77-470a-9512-63fad98dd79a
# ╟─e069c354-8ee3-4490-8937-25a665f0ce54
# ╟─ad99f90c-512f-450e-89e5-3283b6d68a57
# ╟─734e8d32-55de-4eaf-8498-10c9e9f2e51f
# ╟─14c22744-c549-4fac-84e0-7f913f1cfa58
# ╟─ddd0c1bd-19fe-45ab-893b-c94242d612c0
# ╟─a329cfda-d89c-49e0-9f0e-59d2cad465d5
# ╟─88d2172a-9df3-454d-ae8b-290a069a6d27
# ╟─5d1a6abd-3075-4b7f-95a8-2a4101a81185
# ╟─e775138e-2889-4d3b-b9c8-bade656ebc43
# ╟─5232382a-6dd6-42e2-9e9b-92da7dbad672
# ╟─65e570f3-3597-449c-8e8e-7fe563cef9b9
# ╟─ac002d2e-1064-433e-b695-01702083b6a4
# ╟─1cd6bfe2-a720-46a3-adfe-bebfc23189ef
# ╟─816f1e62-3b24-459a-86fb-bfa6624cbe57