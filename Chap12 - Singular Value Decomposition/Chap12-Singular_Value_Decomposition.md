# Chap12

# 특이값 분해(Singular Value Decomposition)



## 12.1 로우-랭크(Low-rank) 행렬에 의한 행렬의 근사

### 12.1.1 로우-랭크 행렬의 이점

Low-rank Assumption은 *'어떤 행렬 $M$ 을 $M$ 보다 더 작은 $rank$ 를 가지는 행렬의 곱으로 표현이 가능하다'* 라는 것을 의미한다.  예를 들어, $rank$ 가 1인 행렬을 생각해 보자. 모든 행들은 1차원 공간에 놓여 있으며, $\{v\}$ 는 이 공간에  대한 기저라고 하자. 행렬의 모든 행은 $v$ 의 스칼라배이다. $u$ 를 이러한 스칼라배들의 벡터라고 하면, 이러한 행렬을  $uv^{T}$ 로 나타낼 수 있다. 이러한 표현으로 작은 저장공간으로 행렬을 나타낼 수 있으므로 로우-랭크 행렬의 이점이라 할 수 있다. 즉, 아래의 식과 같이 랭크 $1$ 인 $m \times n$ 행렬에 대해 단지 $m+n$ 개의 숫자들로 구성할 수 있다.  


$$
\begin{bmatrix} 1 \\ 2 \\ 3 \end{bmatrix}{ \begin{bmatrix} 1 \\ 2 \\ 3 \end{bmatrix} }^{ T }=\begin{bmatrix} 1 \\ 2 \\ 3 \end{bmatrix}\begin{bmatrix} 1 & 2 & 3 \end{bmatrix}=\begin{bmatrix} 1 & 2 & 3 \\ 2 & 4 & 6 \\ 3 & 6 & 9 \end{bmatrix}
$$
마찬가지로, 랭크가 $2$인 행렬은 다음과 같이 표현할 수 있다.


$$
\begin{bmatrix} | & | \\ u_{ 1 } & u_{ 2 } \\ | & | \end{bmatrix}\begin{bmatrix} - & v_{ 1 }^{ T } & - \\ - & v_{ 2 }^{ T } & - \end{bmatrix}
$$
그러나, 실제 데이터에서는 행렬 대부분은 로우-랭크가 아니다. 하지만, 종종 로우-랭크의 근사행렬이 거의 원래 행렬 만큼이나 잘 동작한다고 한다. 

이번 장에서는 주어진 행렬에 대한 최적의 랭크-$k$ 행렬을 찾는 방법에 대해 알아 보도록 한다. 이때의 랭크-$k$ 행렬은 주어진 행렬에 가장 가까운 행렬이다. 



### 12.1.2 행렬의 $Norm$

주어진 행렬에 가장 가까운 랭크-$k$ 행렬을 찾는 문제를 정의하기 위해 행렬들에 대한 거리를 정의하는 것이 필요하다. 벡터들의 경우, 거리는 $norm$ 으로 주어지고, $norm$ 은 내적으로 정의된다. $\mathbb{R}$ 상의 벡터들에 대한 내적은 도트곱이고, 벡터의 $norm$ 은 그 원소들의 제곱의 합의 제곱근이다. 

행렬의 $norm$은 행렬 $A$ 를 벡터로 해석함으로써 계산할 수 있다. $m \times n$ 행렬은 $mn$-벡터로 표현할 수 있다. 

아래와 같이 행렬 $A$ 의 $norm$을 계산하는 것을 *프로베니우스(Frobenius)* $norm$ 이라고 한다.


$$
\left\| A \right\| _{ F }=\sqrt { \sum _{ i }^{  }{ \sum _{ j }^{  }{ A\left[ i,j \right] ^{ 2 } }  }  }
$$
***Lemma*** : $A$의 프로베니우스 $norm$ 의 제곱은 $A$ 의 행들의 제곱의 합과 동일하다.

- **Proof** : $A$ 는 $m \times n$ 행렬이라 하고, $A$ 를 행벡터로 나타내자.  

$$
A=\begin{bmatrix} - & a_{ 1 } & - \\  & \vdots  &  \\ - & a_{ m } & - \end{bmatrix}
$$

- 각각의 행-$i$ 에 대해, 다음과 같이 나타낼 수 있다.

$$
\left\| a_i  \right\| ^{ 2 } = a_i[1]^2 + a_i[2]^2 + \cdots + a_i[n]^2
$$

- 이 식을 프로베니우스 $norm$ 의 식에 대입하면 다음과 같다.

$$
\begin{eqnarray} \left\| A \right\| _{ F } & = & (A[1,1]^{ 2 }+A[1,2]^{ 2 }+\cdots +A[1,n]^{ 2 })+\cdots +(A[m,1]^{ 2 }+A[m,2]^{ 2 }+\cdots +A[m,n]^{ 2 }) \\  & = & \left\| a_{ 1 } \right\| ^{ 2 }+\cdots +\left\| a_{ m } \right\| ^{ 2 } \end{eqnarray}
$$



파이썬의 [`numpy.linalg.norm()`](https://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.linalg.norm.html) 에서는 다양한 $norm$ 계산을 제공하는데 물론 프로베니우스 $norm$ 도 제공한다. 아래의 코드는 `numpy` 모듈을 이용한 프로베니우스 $norm$ 을 구하는 코드이다. 

```python
import numpy as np

A = np.matrix([[1, 0],
               [0, 1]])

fro_norm = np.linalg.norm(A, 'fro')
print(fro_norm)

'''출력결과
1.4142135623730951
'''
```



## 12.2 *트롤리 노선 위치 (Trolley-line-location)* 문제

[9.1 소방차 문제](https://render.githubusercontent.com/view/ipynb?commit=af2d6b76f2034fb89db4a7e8ecc341384cfc48c7&enc_url=68747470733a2f2f7261772e67697468756275736572636f6e74656e742e636f6d2f457863656c73696f72434a482f436f64696e675468654d61747269782f616632643662373666323033346662383964623461376538656363333431333834636663343863372f4368617030392532302d253230546865253230496e6e657225323050726f647563742f4368617030392d5468655f496e6e65725f50726f647563742e6970796e62&nwo=ExcelsiorCJH%2FCodingTheMatrix&path=Chap09+-+The+Inner+Product%2FChap09-The_Inner_Product.ipynb&repository_id=125392345&repository_type=Repository#9.1-%EC%86%8C%EB%B0%A9%EC%B0%A8-%EB%AC%B8%EC%A0%9C) 와 비슷하게 아래의 그림과 같이 벡터 $a_1, ..., a_m$ 으로 표현된 $m$개의 주택의 위치에 대해, 트롤리 노선을 어디에다 위치할 것인지를 찾는 문제를 생각해보자.  

![](./images/trolley.png)



이 문제의 목적은 트롤리 노선을 $m$ 개의 주택에 가능한 가깝게 배치하는 것이다. 

각 벡터 $a_i$ 는 트롤리 노선으로 부터 자신까지의 거리 $d_i$ 를 가진다. 즉, 벡터 $[d_1, ..., d_m]$ 의 $norm$ 을 최소화 해햐한다. 이 문제는 벡터 $norm$ 의 제곱, $d_{1}^{2} + \cdots + d_{m}^{2}$ 을 최소화하는 것과 동일하다.