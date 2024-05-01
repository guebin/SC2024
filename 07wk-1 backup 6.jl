### A Pluto.jl notebook ###
# v0.19.41

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

# ╔═╡ bffeaaa4-f513-4a84-8827-2a69f5c04b29
using Plots

# ╔═╡ f089183d-e83e-4c41-ac17-e48b44126347
md"""
# 07wk-1: 행렬을 보는 관점, 여러 하니
"""

# ╔═╡ 9af868c8-0a3e-4b5b-bae4-7e06f64a7dd6
md"""
## 1. 강의영상
"""

# ╔═╡ 232cfb5c-5921-4759-9632-14d43ff4658b
html"""
<div style="display: flex; justify-content: center;">
<div  notthestyle="position: relative; right: 0; top: 0; z-index: 300;">
<iframe src=
"
https://youtube.com/embed/playlist?list=PLQqh36zP38-w9FQSmInVDk7KkTgK3pPlu&si=bVqOstjkG0MQ1W3B
"
width=600 height=375  frameborder="0" allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
"""

# ╔═╡ 2cbc5887-e09b-411c-bc8f-458831d97f63
md"""
## 2. Imports
"""

# ╔═╡ 718e4953-da21-4136-81dd-b7b3ebf5c038
#using PlutoUI,Images,Plots,Statistics

# ╔═╡ f6239a53-2e36-4aac-b0a2-5af653f5069f
Plots.plotly()

# ╔═╡ 82ce66ae-e37f-4e8d-a316-620985ff798e
PlutoUI.TableOfContents()

# ╔═╡ 2d953e59-1140-440e-8c54-7f0f8a4476e3
md"""
## 3. 생존을 위해 정리한 행렬곱
"""

# ╔═╡ efc41d03-98a0-4de2-af88-c010d8423361
md"""
!!! info "행렬의 벡터화 표현"
	행렬 ${\bf B}_{p \times q}$를 아래와 같이 표한한다고 하자. 
	
	$${\bf B} = \begin{bmatrix} B_1 & B_2 & \dots & B_q \end{bmatrix}$$
	
	여기에서 $B_1,\dots,B_q$는 각각 (p,1)의 차원을 가지는 col-vector 이다. 그러면 행렬 ${\bf A}_{n\times p}$와 행렬 ${\bf B}_{p\times q}$를 곱셉 ${\bf A}{\bf B}$는 아래와 같이 표현가능하다. 
	
	$${\bf A}{\bf B} = \begin{bmatrix}{\bf A}B_1 & {\bf A}B_2 & \dots & {\bf A}B_q\end{bmatrix}$$
"""

# ╔═╡ c4c08183-8ae5-4f3c-9a36-8736b0220123
md"""
-- 예시1: $${\bf B} = \begin{bmatrix} B_1 & B_2 & \dots & B_q \end{bmatrix}$$ 의 표현!
"""

# ╔═╡ 69c67f4b-1ea2-4993-b1b8-76b7aae2bafd
let 
	B1 = [1,2,3]
	B2 = [2,3,4]
	B = [B1 B2]
	B
end 

# ╔═╡ 7eb3c886-2262-419d-83d7-25e2eb1e2e59
md"""
-- 예시2: ${\bf A} {\bf B}$의 표현!
"""

# ╔═╡ 8d12814e-3b25-4768-850c-1bdd9e637e69
let 
	B1 = [1,2,3] # (3,1)
	B2 = [2,3,4] # (3,1)
	B = [B1 B2]	 # (3,2)
	A = [1 2 3; 1 3 4; 1 2 5]
	#A*B
	A * [B1 B2], [A*B1 A*B2]
end 

# ╔═╡ 5b80f068-1dcb-4c88-bbc8-288d136b395d
md"""
## 4. 행렬을 보는 관점
"""

# ╔═╡ b2ae2537-a74d-410a-8bcc-d9f9a5b43b16
md"""
-– 통계학과에서 바라보는 행렬

- 경우1: 데이터의 저장수단 -- 데이터프레임, 이미지, ...
- 경우2: 선형변환 (선형연산)
"""

# ╔═╡ 687c0486-6551-4cb0-9329-2bcad1c99d2f
md"""
### A. 합과 평균
"""

# ╔═╡ 2bd3b211-6306-47ff-a83c-f1c5016a5be9
md"""
-- 아래의 ${\bf X} = \begin{bmatrix} X_1 & X_2 \end{bmatrix}$ 는 데이터를 의미한다.
"""

# ╔═╡ 0ed02da5-b261-4bbc-9549-47f057447d5f
begin
	X1 = randn(100) # (100,1)
	X2 = randn(100) .+ 1 # (100,1)
	X = [X1 X2] # (100,2)
end 

# ╔═╡ e4d88470-25cd-4484-ac32-71a25f03bd87
md"""
-- 아래의 ${\bf j}=[1~ 1~ .... 1]'$는 합계를 구하는 연산을 의미한다.
"""

# ╔═╡ c66e954a-5d68-4c96-9deb-edd79123eca3
j = ones(100) # (100,1)

# ╔═╡ 853b2e64-a587-47c2-8f40-93a41960f65f
j'X1, sum(X1)

# ╔═╡ 8e714a19-9a5f-4749-9b72-b8ddb5f634b6
j'X2, sum(X2)

# ╔═╡ 3918f155-e04c-4bc8-8b34-48ac7f5bb55c
[j'X1 j'X2]

# ╔═╡ 6fa02ff7-9d04-46a4-a513-058c9122287a
j'*[X1 X2]

# ╔═╡ b476c8a9-8c7a-4afe-a6f2-8c7847419258
j'X 

# ╔═╡ 48535f3f-c687-4594-b696-e8b467341e13
md"""
-- $\frac{1}{n}{\bf j}$은 평균연산을 의미한다.
"""

# ╔═╡ 13ae7102-dcbf-4c71-849f-d587db8f47c9
(1/100)j'X1, mean(X1)

# ╔═╡ 2373b3a9-45fa-46ae-acb6-3abfdb48fc51
(1/100)j'X2, mean(X2)

# ╔═╡ 4f06f762-ccdc-499d-a178-ebcf74d517a8
[(1/100)j'X1 (1/100)j'X2]

# ╔═╡ c7b3d727-b767-4665-bf4b-a39af3503a07
(1/100)j'*[X1 X2]

# ╔═╡ ab581c91-69e2-4567-8b7e-ce1acec1bd9d
(1/100)j'X

# ╔═╡ 39b4efd3-88f1-45d1-aee1-93d0173fbd86
md"""
### B. 사영
"""

# ╔═╡ d90ca812-c4c7-420b-bacc-40f9543f95c6
let 
	ϵ = randn(100)
	y = 2X1 .+ 3 + ϵ
	scatter(X1,y,label="(X1,y)")
	X = [j X1]
	H = X*inv(X'X)X' # 변환을 의미하는 행렬, 사영행렬
	scatter!(H*X1, H*y, label="(HX1,Hy)")
	scatter!(X1, H*y, label="(X1,Hy)")	
end 

# ╔═╡ c54b98f7-9216-4cce-baad-824651d819d9
md"""
-- 신기한 거..
"""

# ╔═╡ 1554febf-f082-45ca-9c63-b5c1f55c2e35
let 
	X = [j X1 X2]
	H = X*inv(X'X)X'
	[H*X X]
end 

# ╔═╡ dddbd96c-835a-45f6-9b8c-ff459f48f6a9
md"""
### C. 이미지자료 
"""

# ╔═╡ 28fc34a6-89bc-4337-a8be-e1c78fba9213
md"""
-- RGB 자료형
"""

# ╔═╡ f4f56f7f-6587-406b-8889-cca67d186d1a
RGB(1,0,0)

# ╔═╡ 9a091e7f-0508-4280-a58e-2445d87d6359
md"""
-- 빛의 3원색
"""

# ╔═╡ 1ba51dc3-87d5-4941-a1d0-3925ff6a710b
RGB(1,0,0), RGB(0,1,0), RGB(0,0,1)

# ╔═╡ bef83e0c-210d-45a3-978e-42f014dab06e
RGB(1,1,0), RGB(1,0,1), RGB(0,1,1)

# ╔═╡ 8cef4bb3-e3cd-4a02-807e-70a40f223f4f
RGB(1,1,1), RGB(0,0,0), RGB(0.5,0.5,0.5)

# ╔═╡ a5f0ffc3-9f0c-469e-986e-0cf4b5071704
md"""
-- 색들을 매트릭스로 만들면? 이미지
"""

# ╔═╡ 35d67574-47d4-4001-8052-3f1b0b7bcf8d
[
	RGB(1,0,0) RGB(0,1,0) RGB(0,0,1)
	RGB(1,1,0) RGB(1,0,1) RGB(0,1,1)
	RGB(1,1,1) RGB(0,0,0) RGB(0.5,0.5,0.5)
]

# ╔═╡ c524bd2f-0410-4d59-8613-ff54370c6468
md"""
- 이런식으로 한땀한땀 찍으면 이미지를 만들 수 있다.
"""

# ╔═╡ c82f7cf1-d7be-4ddd-88e4-9fc8983a2644
md"""
-- 프랑스국기
"""

# ╔═╡ a4453c8f-40b4-46f2-b396-f7f3061af2ac
France = [
	RGB(0,0,1) RGB(0,0,1) RGB(1,1,1) RGB(1,1,1) RGB(1,0,0) RGB(1,0,0)
	RGB(0,0,1) RGB(0,0,1) RGB(1,1,1) RGB(1,1,1) RGB(1,0,0) RGB(1,0,0)
	RGB(0,0,1) RGB(0,0,1) RGB(1,1,1) RGB(1,1,1) RGB(1,0,0) RGB(1,0,0)
	RGB(0,0,1) RGB(0,0,1) RGB(1,1,1) RGB(1,1,1) RGB(1,0,0) RGB(1,0,0)
]

# ╔═╡ 32b71fe3-e9ad-437d-907c-035a3a1cb786
md"""
-- 그리스국기
"""

# ╔═╡ d2103a5b-5a16-43dc-9044-648a27ef3bd2
begin
	Greece=load(download("https://upload.wikimedia.org/wikipedia/commons/thumb/5/5c/Flag_of_Greece.svg/1200px-Flag_of_Greece.svg.png?20160309091801"))
	Greece=imresize(Greece, (9,13))
end

# ╔═╡ 024278da-e2ae-4fe8-adbd-20e405efcbdb
Greece[1,1] 

# ╔═╡ 9053aa68-6d3c-4cbf-92d5-09e0d1d781a7
Greece[1:5,1:5] 

# ╔═╡ 0a4d654f-33b5-4b00-a4ce-bc6c1d629ef3
md"""
-- 참고 할 기능: channelview
"""

# ╔═╡ 39353f3c-adac-41b4-a769-07f87f289d67
channelview(France)

# ╔═╡ 8461aa39-cafa-4519-85bc-aee63e1b97a0
channelview(France)[1,:,:] # 빨강

# ╔═╡ 89c0c820-a471-4f2c-b94c-52a49bdc78e7
channelview(France)[2,:,:] # 녹색

# ╔═╡ 753d3cc6-64d3-4f7b-a6dd-f365b17c6b17
channelview(France)[3,:,:] # 파랑

# ╔═╡ 5b6a71c2-bd9c-4bcc-9b93-ef869d425c47
channelview([Greece[1,1]])  # N0f8 무시

# ╔═╡ 2e6c6f2f-5764-4662-aafe-cdf343972c57
RGB(0.051,0.369,0.686), Greece[1,1]

# ╔═╡ 6eaba1ab-be8b-44b0-aaf7-8daea9366296
md"""
-- 참고 할 기능: colorview
"""

# ╔═╡ d5b3df4c-63d0-4f31-82e4-0899516d01fd
colorview(RGB,reshape([0.051,0.369,0.686],3,1,1))

# ╔═╡ 7153503e-ebb0-46a0-9fe8-fd4a62f9082f
colorview(RGB,reshape(rand(3*4*4),3,4,4))

# ╔═╡ 9f709b45-b0c9-436e-af73-5434000ff569
md"""
### D. 컨볼루션
"""

# ╔═╡ df78a889-abb6-4552-9825-07af20bfdec8
md"""
-- 평균을 의미하는 커널
"""

# ╔═╡ 929fc2e7-4391-4d70-a3a4-aaad797b4efe
M = reshape(ones(9)/9,3,3)

# ╔═╡ 72362c9a-51ec-402c-beb5-bdf3c8beb8b0
let 
	X0 = zeros(10,5)
	X1 = ones(10,5)*3
	X = [X0 X1]
	X[3,3] = 27
	X,imfilter(X,M)
	#컨볼루션 = 원소별곱 -> 합-> 이동해서반복
end 

# ╔═╡ fd62e802-c70b-479d-a057-0b5104b5a3b4
France, imfilter(France,M)

# ╔═╡ 5300feeb-ad99-4de2-8585-636f672b1484
Greece, imfilter(Greece,M)

# ╔═╡ 4c906c57-9bb6-4cc3-8fb0-7719cbfed631
md"""
-- 가우시안 커널
"""

# ╔═╡ 115a9087-8eb7-4bd5-a0fe-5b0e3eb8ebc3
G = Kernel.gaussian(1)

# ╔═╡ 17954171-39bf-49f7-8979-d0ab71871002
let 
	X0 = zeros(10,5)
	X1 = ones(10,5)*3
	X = [X0 X1]
	X[3,3] = 27
	X,imfilter(X,G)
	#컨볼루션 = 원소별곱 -> 합-> 이동해서반복
end 

# ╔═╡ 14648835-5bdb-4ef9-b4d2-2af09d7471df
France, imfilter(France,G)

# ╔═╡ ca40ce43-c634-4fb3-b27b-126d63bc4641
Greece, imfilter(Greece,G)

# ╔═╡ 7a25edcb-55d3-4909-b3a8-b73d6ac149a5
md"""
-- 궁금: 왜 가우시안 커널이라고 부를까? 
"""

# ╔═╡ 61f07390-3cbc-4dcd-aee6-fd0bfda51c02
# 가우시안 분포 = 정규분포

# ╔═╡ e956d83c-ed0e-4f15-930e-d4a34d163ae1
surface(Kernel.gaussian(10))

# ╔═╡ c77f4aac-4799-44c2-bda9-7557c95ca143
md"""
## 5. 여러 하니
"""

# ╔═╡ 7a7ace13-924a-4f0f-bcfd-0d6d935a140e
md"""
### A. 하니
"""

# ╔═╡ 2ce4e936-3b78-4729-b52f-cac6f5754849
hani = load(download("https://github.com/guebin/SC2022/blob/main/hani.jpeg?raw=true"))

# ╔═╡ 588dd6c7-177e-4004-8537-6a3664b905c3
size(hani)

# ╔═╡ 3e4d0301-68fc-449d-8400-0779091a11d4
md"""
- 하니는 3024 $\times$ 4032 행렬로 볼 수 있음.
"""

# ╔═╡ 98b84cc5-fa66-414f-9618-3ee6a1daee01
md"""
### B. 하니의 트랜스포즈, 서브매트릭스
"""

# ╔═╡ 8b373d4b-a854-473e-8293-92a002baf567
md"""
-- 하니의 트랜스포즈는 4032 $\times$ 3024 행렬이 된다. 
"""

# ╔═╡ a23ec1b9-f687-4bbd-992d-0355f0330110
hani'

# ╔═╡ 55a1139f-3d9e-4c7b-8507-61fd0a80ff88
md"""
-- 아래는 하니의 서브매트릭스
"""

# ╔═╡ e62d41e7-6006-45e7-ba6e-4018e56f4e65
hani[1000:2000, 1000:2000]

# ╔═╡ 2338452c-64dd-4b35-83be-e0fbba2e090c
hani[1000:1500, 1000:1500]

# ╔═╡ 69763489-0bb1-4beb-92a7-fa84d3a218d5
hani[1200:1300, 1200:1300]

# ╔═╡ 94e8c88a-fda1-468c-9d61-65e2a83d81bd
hani[1200:1205, 1200:1205]

# ╔═╡ 1f97e330-13fa-48f1-b2aa-9cb3a164d9bc
md"""
### C. 루트하니
"""

# ╔═╡ b12fe04c-1bfd-4a0d-872f-79595e808d32
let 
	root_hani = channelview(hani) .|> (x->√x) |> colorview(RGB)
	[hani'  root_hani']
end 

# ╔═╡ b19f2f36-5f7a-4d71-855d-6c130f6e74c1
md"""
-- 값들이 1에 가까워지면서 밝게 보인다. 
- 숫자가 0~1사이이므로 제곱하면 작아지고 루트를 취하면 커진다.
"""

# ╔═╡ 99d83850-0ee4-4a0c-ae90-1870f21c8000
md"""
### D. 후광하니
"""

# ╔═╡ 4f6bc142-590b-4e84-b492-fff84b645095
let 
	shining_hani = channelview(hani) .|> (x-> √x*(x>0.75) + x*(x≤0.75)) |> colorview(RGB)
	[hani' shining_hani']
end 

# ╔═╡ 23c77fb0-770d-4bcd-be19-423afe89ac7c
md"""
### E. 흐릿한 하니
"""

# ╔═╡ d3c1b12d-6f4c-4d24-bde1-76ab491a36ee
md"""
-- 하니에 가우시안커널을 적용하여 스무딩
"""

# ╔═╡ 9dced544-7f75-4da9-a1d7-87f37dfed359
let 
	smoothed_hani = imfilter(hani,Kernel.gaussian(30))
	[hani' smoothed_hani']
end

# ╔═╡ 203493aa-2f0a-4c46-9d4d-8f9063d906e8
md"""
### F. 또렷한 하니
"""

# ╔═╡ 8991ea84-a3d4-427b-bd72-c5d7d64f55df
md"""
-- 개념 
- 이미지 = 스무딩한부분 + 또렷한부분(?) 
- 이미지$\times 2$ =  스무딩한부분$\times 2$ + 또렷한 부분$\times 2$
- 이미지$\times 2$ - 스무딩한부분 = 스무딩한부분 + 또렷한 부분$\times 2$ = 이미지 + 또렷한부분 
"""

# ╔═╡ 43b55707-b168-4684-b84a-2a671f30fd82
md"σ = $(@bind σ Slider(1:50, show_value=true,default=30))"

# ╔═╡ 0394029f-3f5c-453e-8416-40a696e96ecf
let 
	smoothed_hani = imfilter(hani,Kernel.gaussian(σ))
	sharp_hani = hani*2 - smoothed_hani
	[hani' smoothed_hani' sharp_hani']
end

# ╔═╡ 4581a8c3-c858-4503-8e76-6bcd21c4233d
md"""
## 6. 숙제

-- 제곱하니를 만들어보고 원본과 비교해 볼 것
"""

# ╔═╡ Cell order:
# ╟─f089183d-e83e-4c41-ac17-e48b44126347
# ╟─9af868c8-0a3e-4b5b-bae4-7e06f64a7dd6
# ╟─232cfb5c-5921-4759-9632-14d43ff4658b
# ╟─2cbc5887-e09b-411c-bc8f-458831d97f63
# ╠═718e4953-da21-4136-81dd-b7b3ebf5c038
# ╠═bffeaaa4-f513-4a84-8827-2a69f5c04b29
# ╠═f6239a53-2e36-4aac-b0a2-5af653f5069f
# ╠═82ce66ae-e37f-4e8d-a316-620985ff798e
# ╟─2d953e59-1140-440e-8c54-7f0f8a4476e3
# ╟─efc41d03-98a0-4de2-af88-c010d8423361
# ╟─c4c08183-8ae5-4f3c-9a36-8736b0220123
# ╠═69c67f4b-1ea2-4993-b1b8-76b7aae2bafd
# ╟─7eb3c886-2262-419d-83d7-25e2eb1e2e59
# ╠═8d12814e-3b25-4768-850c-1bdd9e637e69
# ╟─5b80f068-1dcb-4c88-bbc8-288d136b395d
# ╟─b2ae2537-a74d-410a-8bcc-d9f9a5b43b16
# ╟─687c0486-6551-4cb0-9329-2bcad1c99d2f
# ╟─2bd3b211-6306-47ff-a83c-f1c5016a5be9
# ╠═0ed02da5-b261-4bbc-9549-47f057447d5f
# ╟─e4d88470-25cd-4484-ac32-71a25f03bd87
# ╠═c66e954a-5d68-4c96-9deb-edd79123eca3
# ╠═853b2e64-a587-47c2-8f40-93a41960f65f
# ╠═8e714a19-9a5f-4749-9b72-b8ddb5f634b6
# ╠═3918f155-e04c-4bc8-8b34-48ac7f5bb55c
# ╠═6fa02ff7-9d04-46a4-a513-058c9122287a
# ╠═b476c8a9-8c7a-4afe-a6f2-8c7847419258
# ╟─48535f3f-c687-4594-b696-e8b467341e13
# ╠═13ae7102-dcbf-4c71-849f-d587db8f47c9
# ╠═2373b3a9-45fa-46ae-acb6-3abfdb48fc51
# ╠═4f06f762-ccdc-499d-a178-ebcf74d517a8
# ╠═c7b3d727-b767-4665-bf4b-a39af3503a07
# ╠═ab581c91-69e2-4567-8b7e-ce1acec1bd9d
# ╟─39b4efd3-88f1-45d1-aee1-93d0173fbd86
# ╠═d90ca812-c4c7-420b-bacc-40f9543f95c6
# ╟─c54b98f7-9216-4cce-baad-824651d819d9
# ╠═1554febf-f082-45ca-9c63-b5c1f55c2e35
# ╟─dddbd96c-835a-45f6-9b8c-ff459f48f6a9
# ╟─28fc34a6-89bc-4337-a8be-e1c78fba9213
# ╠═f4f56f7f-6587-406b-8889-cca67d186d1a
# ╟─9a091e7f-0508-4280-a58e-2445d87d6359
# ╠═1ba51dc3-87d5-4941-a1d0-3925ff6a710b
# ╠═bef83e0c-210d-45a3-978e-42f014dab06e
# ╠═8cef4bb3-e3cd-4a02-807e-70a40f223f4f
# ╟─a5f0ffc3-9f0c-469e-986e-0cf4b5071704
# ╠═35d67574-47d4-4001-8052-3f1b0b7bcf8d
# ╟─c524bd2f-0410-4d59-8613-ff54370c6468
# ╟─c82f7cf1-d7be-4ddd-88e4-9fc8983a2644
# ╠═a4453c8f-40b4-46f2-b396-f7f3061af2ac
# ╟─32b71fe3-e9ad-437d-907c-035a3a1cb786
# ╠═d2103a5b-5a16-43dc-9044-648a27ef3bd2
# ╠═024278da-e2ae-4fe8-adbd-20e405efcbdb
# ╠═9053aa68-6d3c-4cbf-92d5-09e0d1d781a7
# ╟─0a4d654f-33b5-4b00-a4ce-bc6c1d629ef3
# ╠═39353f3c-adac-41b4-a769-07f87f289d67
# ╠═8461aa39-cafa-4519-85bc-aee63e1b97a0
# ╠═89c0c820-a471-4f2c-b94c-52a49bdc78e7
# ╠═753d3cc6-64d3-4f7b-a6dd-f365b17c6b17
# ╠═5b6a71c2-bd9c-4bcc-9b93-ef869d425c47
# ╠═2e6c6f2f-5764-4662-aafe-cdf343972c57
# ╟─6eaba1ab-be8b-44b0-aaf7-8daea9366296
# ╠═d5b3df4c-63d0-4f31-82e4-0899516d01fd
# ╠═7153503e-ebb0-46a0-9fe8-fd4a62f9082f
# ╟─9f709b45-b0c9-436e-af73-5434000ff569
# ╟─df78a889-abb6-4552-9825-07af20bfdec8
# ╠═929fc2e7-4391-4d70-a3a4-aaad797b4efe
# ╠═72362c9a-51ec-402c-beb5-bdf3c8beb8b0
# ╠═fd62e802-c70b-479d-a057-0b5104b5a3b4
# ╠═5300feeb-ad99-4de2-8585-636f672b1484
# ╟─4c906c57-9bb6-4cc3-8fb0-7719cbfed631
# ╠═115a9087-8eb7-4bd5-a0fe-5b0e3eb8ebc3
# ╠═17954171-39bf-49f7-8979-d0ab71871002
# ╠═14648835-5bdb-4ef9-b4d2-2af09d7471df
# ╠═ca40ce43-c634-4fb3-b27b-126d63bc4641
# ╟─7a25edcb-55d3-4909-b3a8-b73d6ac149a5
# ╠═61f07390-3cbc-4dcd-aee6-fd0bfda51c02
# ╠═e956d83c-ed0e-4f15-930e-d4a34d163ae1
# ╟─c77f4aac-4799-44c2-bda9-7557c95ca143
# ╟─7a7ace13-924a-4f0f-bcfd-0d6d935a140e
# ╠═2ce4e936-3b78-4729-b52f-cac6f5754849
# ╠═588dd6c7-177e-4004-8537-6a3664b905c3
# ╟─3e4d0301-68fc-449d-8400-0779091a11d4
# ╟─98b84cc5-fa66-414f-9618-3ee6a1daee01
# ╟─8b373d4b-a854-473e-8293-92a002baf567
# ╠═a23ec1b9-f687-4bbd-992d-0355f0330110
# ╟─55a1139f-3d9e-4c7b-8507-61fd0a80ff88
# ╠═e62d41e7-6006-45e7-ba6e-4018e56f4e65
# ╠═2338452c-64dd-4b35-83be-e0fbba2e090c
# ╠═69763489-0bb1-4beb-92a7-fa84d3a218d5
# ╠═94e8c88a-fda1-468c-9d61-65e2a83d81bd
# ╟─1f97e330-13fa-48f1-b2aa-9cb3a164d9bc
# ╠═b12fe04c-1bfd-4a0d-872f-79595e808d32
# ╟─b19f2f36-5f7a-4d71-855d-6c130f6e74c1
# ╟─99d83850-0ee4-4a0c-ae90-1870f21c8000
# ╠═4f6bc142-590b-4e84-b492-fff84b645095
# ╟─23c77fb0-770d-4bcd-be19-423afe89ac7c
# ╟─d3c1b12d-6f4c-4d24-bde1-76ab491a36ee
# ╠═9dced544-7f75-4da9-a1d7-87f37dfed359
# ╟─203493aa-2f0a-4c46-9d4d-8f9063d906e8
# ╟─8991ea84-a3d4-427b-bd72-c5d7d64f55df
# ╠═43b55707-b168-4684-b84a-2a671f30fd82
# ╠═0394029f-3f5c-453e-8416-40a696e96ecf
# ╟─4581a8c3-c858-4503-8e76-6bcd21c4233d
