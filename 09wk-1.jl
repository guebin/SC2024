### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ 15442a71-00fe-48ed-ae51-def2f2d8876e
using LinearAlgebra, PlutoUI, Plots,Random,CSV,DataFrames

# ╔═╡ e08c8300-db18-11ec-3f60-19c9638fa2d2
md"""
# 09wk-1: 대각화가능행렬
"""

# ╔═╡ 5a141d20-3166-4e3c-96a9-ca5f3395d5f5
md"""
## 1. 강의영상
"""

# ╔═╡ 90db85a9-9864-4e63-98ca-7edb77f1932f
# html"""
# <div style="display: flex; justify-content: center;">
# <div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
# <iframe src=
# "
# https://www.youtube.com/embed/playlist?list=PLQqh36zP38-zjEXJRTGfKk4kWO5T4IlYy
# "
# width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
# """

# ╔═╡ c84fb7f0-e572-46ab-b0a3-7612c1870815
md"""
## 2. Imports
"""

# ╔═╡ 8308a4c4-46b0-4407-aba5-d8b174c2f152
Plots.plotly()

# ╔═╡ 8b5ce365-94e4-4875-ac48-d7d6f4b7409c
PlutoUI.TableOfContents()

# ╔═╡ 5092b326-2d86-410c-8297-baafb75eefb6
md"""
## 3. 귀류법/일반화
"""

# ╔═╡ 59561941-4c46-4380-88ab-b9948b2002ef
md"""
### A. 귀류법
"""

# ╔═╡ 5d1263cf-b4eb-4aa1-a293-e9f8e82fd4b1
md"""
-- 귀류법: 니 논리 대로면… <- 인터넷 댓글에 많음..

```
님 논리대로면..
- XXX가 문제 없으면 서울 전체가 문제가 없고 (애초에 서울은 문제도 아니라는데 왜 이소리는 하고 계신지 모르겠지만)
- 수도권 모 대학이 문제가 없으면 전체가 문제가 없겠네요?
- 지방도 1개 대학이 문제가 없으니 전체가 문제 없겠네요?
와우! 모든 문제가 해결되었습니다! 출산율 감소로 인한 한국대학의 위기가 해결되었.. 아니 애초에 위기가 없었군요!.
어휴.. ㅠㅠ
```
"""

# ╔═╡ f0a0a0af-c54b-4f30-8a5b-f7874ea05e79
md"""
ref: [하이브레인넷](https://m.hibrain.net/braincafe/cafes/38/posts/205/articles/466461?pagekey=466461&listType=TOTAL&pagesize=5&sortType=RDT&limit=30&displayType=TIT&siteid=1&page=1)
"""

# ╔═╡ c85c331a-8cbc-48e8-be68-be0633d4644c
md"""
### B. 일반화
"""

# ╔═╡ c3d6b5bc-bb3e-443a-ab3e-719a9db7bc95
md"""
-- 연필의 정의: 필기도구의 하나. 흑연과 점토의 혼합물을 구워 만든 가느다란 심을 속에 넣고, 겉은 나무로 둘러싸서 만든다. 1565년에 영국에서 처음으로 만들었다.
"""

# ╔═╡ ee58acd1-b1b5-4c90-a392-565bdf778457
md"""
-- 질문: 아래는 연필인가요??
"""

# ╔═╡ 0bac3ee1-5de5-4966-8157-1e85fc5ee7e1
md"""
![](https://guebin.github.io/SP2023/posts/1.%20%EC%B8%A1%EB%8F%84%EB%A1%A0/2023-03-14-2wk-2_files/figure-html/4da5aaa8-0e84-4c7f-a1ae-eae6a2e34076-1-7f6cdfaf-faf5-464e-9ec6-4668e9b1aab1.png)
"""

# ╔═╡ 9d9020be-5c18-48f6-bf87-d3623f568673
md"""
## 4. SVD -- 자투리(minor topics)
"""

# ╔═╡ 05c1b786-159f-4ebd-9748-65c5194133f6
md"""
### A. 비유일성
"""

# ╔═╡ ee1af838-e1ee-4d8d-9dbe-3ec704ea9b43
md"""
!!! info "SVD의 분해결과의 비 유일성" 
	SVD의 분해결과중 ${\bf U}$와 ${\bf V}$는 유일하지 않다. 
"""

# ╔═╡ f3f9396f-7959-4b41-b033-abb0d2514791
md"""
-- 예시1
"""

# ╔═╡ 3007c642-e4bc-4d8c-ba1d-90df5068ee3a
let
	X = [0 1 0 
		 0 0 1 
	     1 0 0]
	U,d,V = svd(X)
	Ũ = -U 
	Ṽ = -V 
	@show Ũ * Ũ' ≈ Ũ' * Ũ ≈ I
	@show Ṽ * Ṽ' ≈ Ṽ' * Ṽ ≈ I
	@show Ũ * Diagonal(d) * Ṽ ≈ X	
end 

# ╔═╡ 02c35de6-acd4-4d5c-ba41-8d7f9523d05a
md"""
-- 예시2
"""

# ╔═╡ 253b6f6a-00c2-4b4f-a2b5-b651a9b957ef
let
	X = [0 1 0 
		 0 0 1 
	     1 0 0]
	U,d,V = svd(X)
	U1,U2,U3 = eachcol(U)
	V1,V2,V3 = eachcol(V)
	Ũ = [U1 U2 -U3] 
	Ṽ = [V1 V2 -V3] 
	@show Ũ * Ũ' ≈ Ũ' * Ũ ≈ I
	@show Ṽ * Ṽ' ≈ Ṽ' * Ṽ ≈ I
	@show Ũ * Diagonal(d) * Ṽ ≈ X
end 

# ╔═╡ 3435be8f-5760-495b-8389-94516b37ab98
md"""
### B. $d_i$는 음수가 아님.
"""

# ╔═╡ 6a141456-80f8-4e0a-a229-90a79f6268de
md"""
!!! info "SVD에서 D의 대각원소는 음수가 아님" 
	SVD의 분해결과중 ${\bf D}$의 대각원소는 음수가 아니도록 "정의"되어있다. 원소가 0이 되는 경우는 있을 수 있다. 
"""

# ╔═╡ 0572eeb0-07c2-45f0-afd4-dae608713903
md"""
-- 예시1
"""

# ╔═╡ 582c59f4-fa61-4a47-bbc5-44056af90f81
let
	X = [0 1 0 
		 0 0 1 
	     1 0 0]
	U,d,V = svd(X)
	@show d
end 

# ╔═╡ d24a006a-7c28-4d3b-a89d-7d5d0eb4cacd
md"""
-- 예시2
"""

# ╔═╡ bb3ca1cf-e053-422b-9eab-691d3895b5cf
let
	X = [2  0
		 0 -1]
	λ,Ψ = eigen(X)
	U = V = Ψ
	d = λ; D=Diagonal(d)
	@show d
	@show U * U' ≈ U' * U ≈ I
	@show V * V' ≈ V' * V ≈ I # 이건 사실 안해봐도 되는데,.	
	@show U * D * V' ≈ X 
	println("--여기에서의 U,D,V 는 X=UDV' 를 만족하지만 SVD의 분해결과로 인정하지 않음--\n")
	U1,U2 = eachcol(U)
	d1,d2 = d
	Ũ = [-U1 U2]
	d̃ = [-d1, d2]
	D̃ = Diagonal(d̃)
	@show d̃
	@show Ũ * Ũ' ≈ Ũ' * Ũ ≈ I
	@show Ũ * D̃ * V' ≈ X 
	println("--여기에서의 Ũ,D̃,V 는 X=UDV' 를 만족하지만 SVD의 분해결과로 인정--")
end 

# ╔═╡ 154b19b3-fcf4-4db4-a36f-9c862c4b386e
md"""
### C. 고유분해와 특이값분해의 관계
"""

# ╔═╡ e65f9768-a3b9-4f6d-8751-0829205b0089
md"""
-- 경험1: 고유분해와 특이값분해의 결과가 일치하는 경우가 있었음.
"""

# ╔═╡ b74f600b-001c-4749-a93d-c9d28362a723
md"""
-- 경험2: 주성분분석을 특이값분해/고유뷴해 양쪽을 이용해서 계산할 수 있었음. 
"""

# ╔═╡ 1ed0e5d5-bc1a-4594-82eb-f6311cbbb663
md"""
-- 본질: 사실 특이값분해는 고유분해에서 유도된 것임~!
"""

# ╔═╡ 05ee081c-a361-4d3e-94df-12d9af411779
md"""
앞으로 자세히 알아봅시다
"""

# ╔═╡ 006740e0-f12c-4312-b181-30d5e6f82ec3
md"""
## 5. 고유값과 고유벡터
"""

# ╔═╡ cba26403-d491-4eae-be95-7a11984c2b0e
md"""
### A. 정의
"""

# ╔═╡ 82dc78d3-4dea-47c4-bf93-a5e162a89916
md"""
!!! info "고유값과 고유벡터의 정의"
	임의의 정사각행렬 ${\bf A}_{n\times n}$에 대하여 어떠한 벡터 ${\psi}_{n\times 1}\neq 0$ 가 적당한한 값 $\lambda$에 대하여 

	${\bf A}{\psi} = \lambda \psi$

	를 만족하면 $\psi$를 $\lambda$의 고유벡터라고 하고 $\lambda$는 $\psi$에 대응하는 고유값이라고 한다. 
"""

# ╔═╡ 7b3b6699-5666-4262-a648-f62cb20627ad
md"""
!!! info "고유값과 고유벡터의 정의"
	0-벡터는 고유벡터로 인정하지 않는다는 것을 주의할 것!!
"""

# ╔═╡ 1e655ddb-83bf-4b41-b0d0-e57f957cfc25
md"""
### B. 구하는 방법 1
"""

# ╔═╡ fe4085be-e310-4b92-889d-7e5fd8006bbb
md"""
-- 그냥 찾음.. 
"""

# ╔═╡ af49f4f4-86a2-40dd-970e-df7a0183e5e3
md"""
-- 예제1: 단위행렬의 고유값과 고유벡터를 찾아라.
"""

# ╔═╡ bbae5460-6afd-4d46-99dc-6f8c283e7ffb
md"""
-- 예제2: 0행렬의 고유값과 고유벡터를 찾아라.
"""

# ╔═╡ f8044f4f-9269-4955-ade6-0ec5f6889412
md"""
-- 예제3: 대각행렬의 고유값과 고유벡터를 찾아라.
"""

# ╔═╡ 62015917-017d-457d-a975-2eff5f581b3b
md"""
### C. 구하는 방법 2 -- 특성방정식
"""

# ╔═╡ fe5876ba-30f9-470d-bd42-4c188fa580f3
md"""
!!! info "고유값/고유벡터를 체계적으로 찾는방법"
	임의의 정사각행렬 ${\bf A}$의 고유값은 항상 아래를 만족하는 $\lambda$를 풀어서 찾을 수 있다. 
	
	$\det({\bf A}-\lambda {\bf I})=0$

	하나의 고유값이 정해지면 그 이후에는 아래의 식을 풀어서 $\psi$를 찾을 수 있다. 

	${\bf A}{\psi} = \lambda \psi$
"""

# ╔═╡ bd9fd3d7-0303-4098-afb6-c169e36ed089
md"""
### D. 고유값/고유벡터에 대한 고찰
"""

# ╔═╡ 8aeab63a-6c60-4d4b-9083-b89e96232973
md"""
!!! warning "고유값은 항상 존재 + 고유값은 중복될 수 있음 + 고유값이 실수라는 보장없음"
	임의의 $2\times 2$ 행렬 ${\bf A}$ 를 고려하자. 

	$\det({\bf A}-\lambda {\bf I})=0$

	는 $\lambda$에 대한 이차방정식의 형태일 것이다. 아래의 사실을 쉽게 유추할 수 있다. 

	1. 이차방정식의 근은 항상 존재하므로 ${\bf A}$의 고유값은 항상 존재한다. (그리고 고유값/고유벡터의 정의에 따라서 ${\bf A}$에 대응하는 고유벡터도 항상 존재한다) 
	2. 이차방정식은 중근을 가질 수 있으므로 ${\bf A}$가 서로다른 두개의 고유값을 가진다는 보장은 없다. 
	3. 또한 이차방정식은 허근을 가질 수 있으므로 ${\bf A}$의 고유값이 실수라는 보장 역시 없다. 

	이러한 논의는 임의의 $n \times n$ 행렬 ${\bf A}$ 에 대하여도 확장된다. 즉 아래가 성립한다. 

	1. ``{\bf A}``의 고유값과 고유벡터는 항상 존재한다. (대수학의 기본정리) --> 좋은데?
	2. ``{\bf A}``가 서로 다른 $n$개의 고유값을 가진다는 보장은 없다. --> rank???
	3. ``{\bf A}``의 고유값이 실수라는 보장은 없다. --> 허수같은건 생각하기 싫은걸??
"""

# ╔═╡ ca93f2b5-0117-48cd-87f0-85f3f9d2ca23
md"""
!!! warning "고유벡터의 비 유일성"
	고유벡터는 언제나 유일하지 않다. 예를들어 $\psi$가 고유벡터의 정의를 만족한다면 $-\psi, \sqrt{2}\psi$ 역시 고유벡터의 정의를 만족한다. 즉 고유벡터의 정의를 만족하면 $\psi$를 하나만 찾으면 정의를 만족하는 ${\boldsymbol \psi}$ 무한개를 찾을 수 있다. --> 싫은걸? 표준화 하고 싶은걸..?
"""

# ╔═╡ 68ec2ffe-cfd7-43ed-b2cc-23ddf0245aa1
md"""
!!! info "고유값과 고유벡터의 정의 -- 매트릭스 버전"
	임의의 정사각행렬 ${\bf A}_{n\times n}$에 대한 고유값들을 $\lambda_1,\dots,\lambda_n$ 이라고 하고 그것에 대응하는 고유벡터들을 ${\boldsymbol \psi}_1\dots {\boldsymbol \psi}_n$ 이라고 하자. 그리고 
	- ``{\boldsymbol \Psi} = \begin{bmatrix} {\boldsymbol \psi}_1 & {\boldsymbol \psi}_2 & \dots& {\boldsymbol \psi}_n \end{bmatrix}``
	- ``{\boldsymbol \Lambda} = \text{diag}(\lambda_1, \dots,\lambda_n)``
	와 같은 매트릭스를 정의하자. 그렇다면 아래의 식이 만족한다. 

	${\bf A}{\boldsymbol \Psi} = {\boldsymbol \Psi}{\boldsymbol \Lambda}$

	이러한 표현은 임의의 정사각행렬 ${\bf A}$에 대하여 언제나 가능함을 기억하자. 그리고 앞으로는 편의를 위하여 ``{\boldsymbol \Lambda}``를 고유값 행렬이라고 하고 ${\boldsymbol \Psi}$ 를 고유벡터행렬이라고 부르자.
"""

# ╔═╡ 9f3e1a8f-95d1-47ca-bf67-9308bb224c10
md"""
### E. 쓸모 없는 내용..
"""

# ╔═╡ 503c9c66-b1a9-4bf3-813a-980635491ce2
md"""
-- 원래 고유값이 존재하면 그에 대응하는 고유벡터도 당연히 존재한다. (고유값과 고유벡터는 홀로 정의되는게 아님, pair로 정의됨)
"""

# ╔═╡ 3634c326-353f-43c8-8170-2efea684a8cb
md"""
-- 그런데 "특성방정식의 근이 존재한다면, 대응하는 고유벡터가 존재한다" 라는게 직관적으로 와닿지 않는 경우가 있음. 이를 위해서 "특성방정식의 근에 대응하는 고유벡터가 항상 존재함"을 증명해보자. 
"""

# ╔═╡ 37df0670-5120-4eb6-86bc-31f7d0864fef
md"""
!!! info "선형독립 -- 벡터버전"
	벡터 $A_1,A_2,\dots,A_n$를 고려하자. 만약 

	$\psi_1 A_1 + \psi_2A_2 + \dots + \psi_nA_n=0$

	을 만족하는 경우가 $\psi_1 = \psi_2 = \dots = \psi_n=0$ 뿐이라면 $A_1,A_2,\dots,A_n$ 은 선형독립이고, 그렇지 않다면 $A_1,A_2,\dots,A_n$ 은 선형독립이 아니다.
	
"""

# ╔═╡ b2b8fb13-2427-4f27-9f1e-76321504698b
md"""
!!! info "선형독립 -- 매트릭스 버전"
	벡터 $A_1,A_2,\dots,A_n$을 column으로 가지는 매트릭스 ${\bf A}=[A_1 ~ A_2~ \dots ~ A_n]$ 를 고려하자. 만약 
	
	$${\bf A}{\boldsymbol \psi}=0$$ 

	만족하는 ${\boldsymbol \psi}$가 0-벡터 밖에 없을 경우는 ${\bf A}$의 column들이 선형독립이고 (즉 벡터 $A_1,\dots,A_n$은 선형독립이고), 그렇지 않다면 (= 위의식을 만족하는 ${\boldsymbol \psi}$가 0-벡터 이외에 존재하면) ${\bf A}$의 column들은 선형독립이 아니다.
"""

# ╔═╡ d6de6890-5798-4ccf-af44-4af4edf4ba7e
md"""
!!! caution "특성방정식의 근에 대응하는 고유벡터"
	특성방정식의 근 $\lambda$에 대응하는 고유벡터가 반드시 하나는 존재한다. 
	
	(proof) 귀류법을 쓰자. 즉 "고정된 $\lambda^*$에 대응하는 고유벡터가 없다"고 가정하자. 편의상 특성방정식을 만족하는 하나의 근 $\lambda^*$를 fix할 때, $\lambda^*$에 대응하는 고유벡터가 없다는 의미는 

	$$({\bf A}-\lambda^*{\bf I}){\boldsymbol \psi} =0$$

	를 만족하는 ${\boldsymbol \psi}$는 오직 ${\boldsymbol \psi}={\bf 0}$ 뿐이라는 것을 의미한다. 즉 $({\bf A}-\lambda^*{\bf I})$ 는 선형독립이라는 의미이다. 그런데 이는 사실이 아니다. 왜냐하면 $\lambda^*$는 

	$$\det({\bf A}-\lambda^*{\bf I})=0$$

	를 만족하고 따라서 행렬 ${\bf A}-\lambda^*{\bf I}$는 역행렬이 없는 행렬이 되고 이렇게 되면 ${\bf A}-\lambda^*{\bf I}$의 column 들은 선형독립이 아니게 된다. (모순)
"""

# ╔═╡ ef80f82f-6acf-4848-92d0-56f554670055
md"""
## 6. 대각화가능행렬
"""

# ╔═╡ 2b4c2b61-53b6-4585-9e5e-774a0cb9bc9f
md"""
### A. 고유벡터의 랭크
"""

# ╔═╡ 9d71c79e-f1f5-4253-881e-44e6d89fef8a
md"""
-- 고유벡터가 유일하지 않다는 것은 이미 알고있는 사실임.
"""

# ╔═╡ 23729642-3a1e-4f22-b05b-356fad1bf54a
md"""
-- 예제1: 소스(?)가 2개인 경우
"""

# ╔═╡ bc325027-c401-410d-be4e-2e99e5077269
let 
	A = [1 0 
		 0 2]
	eigen(A)
end

# ╔═╡ 73808ba1-63b7-4157-b606-26814a864f2f
md"""
-  $\begin{bmatrix} 1 \\ 0 \end{bmatrix}$ 이 고유벡터 이므로 $\begin{bmatrix} 0.5 \\ 0 \end{bmatrix}$, $\begin{bmatrix} -3.14 \\ 0 \end{bmatrix}$ 등도 고유벡터이다. 
"""

# ╔═╡ 2243a0ee-3329-4b06-883c-6dd8bfd44551
md"""
-  $\begin{bmatrix} 0 \\ 1 \end{bmatrix}$ 이 고유벡터 이므로 $\begin{bmatrix} 0 \\ 5.3 \end{bmatrix}$, $\begin{bmatrix} 0 \\ -\sqrt{3} \end{bmatrix}$ 등도 고유벡터이다. 
"""

# ╔═╡ 14c834bd-61ef-45d8-8aa1-cfee057f1e3d
md"""
-- 예제2: 소스(?)가 1개인 경우
"""

# ╔═╡ cbdabff5-c038-4f7a-ac35-545aa7127bd3
let 
	B = [1 1 
		 0 1]
	eigen(B)
end

# ╔═╡ 212b9778-dee0-45a2-a4a8-b75f4fa7a0bb
md"""
-  $\begin{bmatrix} 1 \\ 0 \end{bmatrix}$ 이 고유벡터 이므로 $\begin{bmatrix} -1 \\ 0 \end{bmatrix}$, $\begin{bmatrix} -3.14 \\ 0 \end{bmatrix}$ 등도 고유벡터이다. 
"""

# ╔═╡ 8657dff6-92e6-45df-853d-d13ee2959138
md"""
-- 예제1,2 모두 무한개의 고유벡터를 가지는 것은 맞지만 느낌이 다르다. 

-- 예제2의 경우 $\begin{bmatrix} 1 \\ 0\end{bmatrix}$을 소스로 한 벡터들의 조합을 만들 수 있지만 예제1의 경우 $\begin{bmatrix} 1 \\ 0\end{bmatrix}$, $\begin{bmatrix} 0 \\ 1\end{bmatrix}$를 소스로 한 조합을 만들 수 있기 때문이다.

-- 즉 예제1이 가질 수 있는 고유벡터 조합들이 예제2보다 더 풍부하다. 이 풍부한 느낌을 좀 더 수학적으로 표현할 수 있을까? $\to$ 이게 바로 rank의 개념이죠??
"""

# ╔═╡ 9f690cb6-a65f-41c3-a3f1-fb2000d54500
md"""
### B. Full-rank 고유벡터 행렬
"""

# ╔═╡ 4ab074e8-af46-4db4-9cf3-23a4bfa335a9
md"""
!!! info "대각화가능 (고유벡터행렬의 full-rank matrix 일때)"
	임의의 정사각행렬 ${\bf A}_{n\times n}$에 대한 아래를 만족하는 고유값행렬과 고유벡터행렬은 항상 존재한다. 
	
	${\bf A}{\boldsymbol \Psi} = {\boldsymbol \Psi}{\boldsymbol \Lambda}$

	만약에 ${\boldsymbol \Psi}$가 full-rank-matrix 라면 (= ${\boldsymbol \Psi}^{-1}$ 이 존재한다면) 아래와 같은 표현이 가능하다. 

	${\bf A} = {\boldsymbol \Psi}{\boldsymbol \Lambda}{\boldsymbol \Psi}^{-1}$

	그리고 이때 ${\bf A}$는 **"대각화가능"**행렬이라고 부른다. (수식 ${\bf \Lambda} = {\boldsymbol \Psi}^{-1}{\bf A}{\boldsymbol \Psi}$ 를 관찰해보시면 대각화가능행렬이라는 이름이 생긴 이유를 알 수 있습니다!~) 
"""

# ╔═╡ c0b4b624-23d4-4d6d-a69f-df6e8a77b69a
md""" 
!!! warning "대각화가능행렬 = 고유분해 표현이 존재하는 행렬"
	아래의 수식
	
	${\bf A} = {\boldsymbol \Psi}{\boldsymbol \Lambda}{\boldsymbol \Psi}^{-1}$
	
	를 ${\bf A}$의 고유분해 (eigendecomposition) 라고 표현한다. 따라서 아래의 용어들은 문맥상 같은 말이다. 

	- ``{\bf A}`` 가 대각화가능하다. 
	- ``{\bf A}`` 의 고유분해표현이 존재한다. 
	- ``{\bf A}`` 를 고유분해할 수 있다. 

	통계학과에서는 고유분해라는 표현을 자주쓰므로 사실 ``{\bf A}``를 "고유분해 가능행렬" 정도의 용어로 정의하고 싶은데 아쉽게도 그러한 표현은 쓰지 않는다. (쓰는거 보면 알려주세요..)
"""

# ╔═╡ f6bf3f99-485a-459b-94f6-416b7c188f85
md"""
### C. 대각화가능 여부 체크하기
"""

# ╔═╡ 278afdaf-26f1-44d1-9954-a54dc498c257
md"""
!!! info "Fact1"
	서로 다른 $k$개의 고유값에 대응하는 고유벡터행렬은 rank $k$ 이다. 즉 full-rank 이다.  
"""

# ╔═╡ ab4bed72-dd05-4279-9383-92a51a61d17b
md"""
!!! info "Fact2"
	고유값이 중근을 가진다면 (근이 2개 겹친다면), 대응하는 고유벡터행렬의 rank는 $1$일수도 있고 $2$일수도 있다. 만약에 고유값이 $m$개의 중첩된 근을 가진다면 고유벡터행렬의 rank는 최소 $1$ 최대 $m$이다. 
"""

# ╔═╡ b93d4e26-a6d7-4a56-8c12-ef57eef6ba3e
md"""
!!! info "대각화가능여부를 체크하는 팁? -- 안중요한데 모르면 좀 아쉬울때 있음"

	고유벡터행렬을 구한뒤에 rank를 조사하는 것이 가장 깔끔하지만, 그게 어려울 때도 있다. 

	1. 모든 고유값이 서로 다르다면 대각화가능하다. 
	2. 중복된 고유값이 있다면, 그 중복된 고유값들에 대응하는 고유벡터행렬이 full-rank 인지만 조사하면 대각화가능하다고 판단할 수 있다. 
"""

# ╔═╡ a2f29756-d44c-4536-a505-8752fe8e805a
md"""
### D. 다음시간 예고편: 넘어야할 산..
"""

# ╔═╡ 36a7faf7-3d85-47ed-9f40-b88d8364abc8
md"""
!!! info "대각화가능 -> ... -> 고유값분해 = 특이값분해"
	임의의 정사각행렬 ${\bf A}_{n\times n}$가 대각화가능행렬이라고 하자. 즉 아래를 만족하는 고유값행렬과 고유벡터행렬이 존재한다. 

	${\bf A} = {\boldsymbol \Psi}{\boldsymbol \Lambda}{\boldsymbol \Psi}^{-1}$

	여기에서 **우연히** 

	1. ``{\bf \Psi}^\top= {\bf \Psi}^{-1}`` 이 성립하고 (즉 ${\bf \Psi}$ 가 직교행렬이고)
	2. 모든 고유값이 양수라면

	아래와 같이 표현가능한 직교행렬 ${\boldsymbol \Psi}={\bf U}={\bf V}$ 와 대각행렬 ${\bf D}={\bf \Lambda}$를 항상 정의할 수 있다. 

	${\bf A} = {\boldsymbol \Psi}{\boldsymbol \Lambda}{\boldsymbol \Psi}^\top = {\bf U} {\bf D} {\bf V}^\top$

	즉 이 경우 고유값분해=특이값분해가 된다. 
	
"""

# ╔═╡ de97e015-f2d8-453e-9f71-e748eab4159e
md"""
## 7. 숙제 

정사각행렬 ${\bf A}_{2\times 2}$가 0행렬인 경우에 ${\bf A}$의 고유벡터들이 span하는 공간이 몇차원인가? 
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[compat]
CSV = "~0.10.14"
DataFrames = "~1.6.1"
Plots = "~1.40.4"
PlutoUI = "~0.7.58"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.2"
manifest_format = "2.0"
project_hash = "e32772c00eac22223105502786c1f8baa8386ce2"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "297b6b41b66ac7cbbebb4a740844310db9fd7b8c"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.1"

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

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "6c834533dc1fabd820c1db03c839bf97e45a3fab"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.14"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a4c43f59baa34011e303e76f5c8c91bf58415aaf"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+1"

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
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

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
git-tree-sha1 = "6cbbd4d241d7e6579ab354737f4dd95ca43946e1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.1"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "97d79461925cdb635ee32116978fc735b9463a39"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.19"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

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

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "9f00e42f8d99fdde64d40c8ea5d14269a2e2c1aa"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.21"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

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
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

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

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "3437ade7073682993e092ca570ad68a2aba26983"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.3"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a96d5c713e6aa28c242b0d25c1347e258d6541ab"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.3+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "359a1ba2e320790ddbe4ee8b4d54a305c0ea2aff"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.0+0"

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
git-tree-sha1 = "8e59b47b9dc525b70550ca082ce85bcd7f5477cd"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.5"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

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

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

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
git-tree-sha1 = "e0b5cd21dc1b44ec6e64f351976f961e6f31d6c4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.3"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dae976433497a2f841baadea93d27e68f1a12a97"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.39.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0a04a1318df1bf510beb2562cf90fb0c386f58c4"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.39.3+1"

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
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

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
git-tree-sha1 = "3da7367955dcc5c54c1ba4d402ccdc09a1a3e046"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+1"

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
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "442e1e7ac27dd5ff8825c3fa62fbd1e86397974b"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.4"

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

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

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

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "88b895d13d53b5577fd53379d913b9ab9ac82660"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "0e7508ff27ba32f26cd459474ca2ede1bc10991f"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.1"

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
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

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
git-tree-sha1 = "71509f04d045ec714c4748c785a59045c3736349"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.7"
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

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "532e22cf7be8462035d092ff21fada7527e2c488"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.6+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

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
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

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
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

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
# ╟─e08c8300-db18-11ec-3f60-19c9638fa2d2
# ╟─5a141d20-3166-4e3c-96a9-ca5f3395d5f5
# ╠═90db85a9-9864-4e63-98ca-7edb77f1932f
# ╟─c84fb7f0-e572-46ab-b0a3-7612c1870815
# ╠═15442a71-00fe-48ed-ae51-def2f2d8876e
# ╠═8308a4c4-46b0-4407-aba5-d8b174c2f152
# ╠═8b5ce365-94e4-4875-ac48-d7d6f4b7409c
# ╟─5092b326-2d86-410c-8297-baafb75eefb6
# ╟─59561941-4c46-4380-88ab-b9948b2002ef
# ╟─5d1263cf-b4eb-4aa1-a293-e9f8e82fd4b1
# ╟─f0a0a0af-c54b-4f30-8a5b-f7874ea05e79
# ╟─c85c331a-8cbc-48e8-be68-be0633d4644c
# ╟─c3d6b5bc-bb3e-443a-ab3e-719a9db7bc95
# ╟─ee58acd1-b1b5-4c90-a392-565bdf778457
# ╟─0bac3ee1-5de5-4966-8157-1e85fc5ee7e1
# ╟─9d9020be-5c18-48f6-bf87-d3623f568673
# ╟─05c1b786-159f-4ebd-9748-65c5194133f6
# ╟─ee1af838-e1ee-4d8d-9dbe-3ec704ea9b43
# ╟─f3f9396f-7959-4b41-b033-abb0d2514791
# ╠═3007c642-e4bc-4d8c-ba1d-90df5068ee3a
# ╟─02c35de6-acd4-4d5c-ba41-8d7f9523d05a
# ╠═253b6f6a-00c2-4b4f-a2b5-b651a9b957ef
# ╟─3435be8f-5760-495b-8389-94516b37ab98
# ╟─6a141456-80f8-4e0a-a229-90a79f6268de
# ╟─0572eeb0-07c2-45f0-afd4-dae608713903
# ╠═582c59f4-fa61-4a47-bbc5-44056af90f81
# ╟─d24a006a-7c28-4d3b-a89d-7d5d0eb4cacd
# ╠═bb3ca1cf-e053-422b-9eab-691d3895b5cf
# ╟─154b19b3-fcf4-4db4-a36f-9c862c4b386e
# ╟─e65f9768-a3b9-4f6d-8751-0829205b0089
# ╟─b74f600b-001c-4749-a93d-c9d28362a723
# ╟─1ed0e5d5-bc1a-4594-82eb-f6311cbbb663
# ╟─05ee081c-a361-4d3e-94df-12d9af411779
# ╟─006740e0-f12c-4312-b181-30d5e6f82ec3
# ╟─cba26403-d491-4eae-be95-7a11984c2b0e
# ╟─82dc78d3-4dea-47c4-bf93-a5e162a89916
# ╟─7b3b6699-5666-4262-a648-f62cb20627ad
# ╟─1e655ddb-83bf-4b41-b0d0-e57f957cfc25
# ╟─fe4085be-e310-4b92-889d-7e5fd8006bbb
# ╟─af49f4f4-86a2-40dd-970e-df7a0183e5e3
# ╟─bbae5460-6afd-4d46-99dc-6f8c283e7ffb
# ╟─f8044f4f-9269-4955-ade6-0ec5f6889412
# ╟─62015917-017d-457d-a975-2eff5f581b3b
# ╟─fe5876ba-30f9-470d-bd42-4c188fa580f3
# ╟─bd9fd3d7-0303-4098-afb6-c169e36ed089
# ╟─8aeab63a-6c60-4d4b-9083-b89e96232973
# ╟─ca93f2b5-0117-48cd-87f0-85f3f9d2ca23
# ╟─68ec2ffe-cfd7-43ed-b2cc-23ddf0245aa1
# ╟─9f3e1a8f-95d1-47ca-bf67-9308bb224c10
# ╟─503c9c66-b1a9-4bf3-813a-980635491ce2
# ╟─3634c326-353f-43c8-8170-2efea684a8cb
# ╟─37df0670-5120-4eb6-86bc-31f7d0864fef
# ╟─b2b8fb13-2427-4f27-9f1e-76321504698b
# ╟─d6de6890-5798-4ccf-af44-4af4edf4ba7e
# ╟─ef80f82f-6acf-4848-92d0-56f554670055
# ╟─2b4c2b61-53b6-4585-9e5e-774a0cb9bc9f
# ╟─9d71c79e-f1f5-4253-881e-44e6d89fef8a
# ╟─23729642-3a1e-4f22-b05b-356fad1bf54a
# ╠═bc325027-c401-410d-be4e-2e99e5077269
# ╟─73808ba1-63b7-4157-b606-26814a864f2f
# ╟─2243a0ee-3329-4b06-883c-6dd8bfd44551
# ╟─14c834bd-61ef-45d8-8aa1-cfee057f1e3d
# ╠═cbdabff5-c038-4f7a-ac35-545aa7127bd3
# ╟─212b9778-dee0-45a2-a4a8-b75f4fa7a0bb
# ╟─8657dff6-92e6-45df-853d-d13ee2959138
# ╟─9f690cb6-a65f-41c3-a3f1-fb2000d54500
# ╟─4ab074e8-af46-4db4-9cf3-23a4bfa335a9
# ╟─c0b4b624-23d4-4d6d-a69f-df6e8a77b69a
# ╟─f6bf3f99-485a-459b-94f6-416b7c188f85
# ╟─278afdaf-26f1-44d1-9954-a54dc498c257
# ╟─ab4bed72-dd05-4279-9383-92a51a61d17b
# ╟─b93d4e26-a6d7-4a56-8c12-ef57eef6ba3e
# ╟─a2f29756-d44c-4536-a505-8752fe8e805a
# ╟─36a7faf7-3d85-47ed-9f40-b88d8364abc8
# ╟─de97e015-f2d8-453e-9f71-e748eab4159e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
