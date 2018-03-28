# Chap13 

# 고유벡터 - EigenVector



*13.1 ~ 13.2 생략*



## 13.3 고유값과 고유벡터

***Definition*** : 정방행렬 $A$에 대하여, 스칼라(scalar)인 $\lambda$와 영이 아닌 벡터 $v$에 대해 $Av=\lambda v$가 만족하는 경우, $\lambda$는 $A$의 *고유값(eigenvalue)*,  $v$ 는 대응하는 *고유벡터(eigenvector)* 라고 한다.



만약, $\lambda$가 행렬 $A$의 고유값이면, 대응하는 고유벡터는 무수히 많다. 집합 $\{v : Av = \lambda v \}$ 는 벡터공간이며 고유값 $\lambda$에 대응하는 *고유공간(eigenspace)* 이라 한다. 따라서, 고유공간에 있는 임의의 영이 아닌 벡터는 고유벡터로 간주된다. 일반적으로 고유벡터의 크기($norm$)가 1이라는 제한을 두는 것이 다루기에 편리하다.



***Example 13.3.3*** : 행렬 $A$를 아래와 같은 대각행렬이라 하자.
$$
A=\begin{bmatrix} \lambda _{1} &  &  \\  & \ddots  &  \\  &  & \lambda _{ n } \end{bmatrix}
$$
행렬 $A$에 대한 고유벡터와 고유값은 무엇일까? 표준 기저 벡터 $e_1, \dots, e_n$ 에 대해, $Ae_1 = \lambda_1 , \dots, Ae_n = \lambda_n e_n$ 이므로, $e_1, \dots, e_n$ 은 고유벡터이고 대각원소인 $\lambda_1, \dots , \lambda_n$ 은 고유값이다.



***Example 13.3.5*** : 행렬 $A$의 한 고유값은 0이라고 하자. 이 고유값에 대응하는 고유벡터는 $Av = 0v$를 만족하는 영이 아닌 벡터 $v$ 이다. 즉, 벡터 $v$ 는 $Av$가 영벡터가 되게 하는 영이 아닌 벡터이며, $v$ 는 영공간(null space)에 속한다. 역으로, 만약 $A$의 영공간이 자명하지 않으면 $0$은 $A$의 고유값이다.

위의 Example 13.3.5는 고유값 0에 대응하는 고유벡터(즉, 영공간에 속하는 영이 아닌 벡터)를 찾는 방법에 대한 설명이다.  행렬 $A$ 의 고유값을 $\lambda$, 대응하는 고유벡터를 $v$라고 하면, $Av = \lambda v$ 이다. 따라서, $Av - \lambda v = 0$ 이다. $Av - \lambda v = (A - \lambda I)v $ 이므로  $(A - \lambda I)v$ 는 영벡터이다. 이것은 벡터 $v$가 $A - \lambda I$의 영공간에 속하는 영이 아닌 벡터임을 의미한다. 따라서, $A - \lambda I$는 비가역적이다.



***Corollary*** : 만약 $\lambda$ 가 행렬 $A$의 고유값일 경우 $\lambda$ 는 또한 $A^T$의 고유값이다.

- **Proof** : $\lambda$를 행렬 $A$의 고유값이라 하면, $A-\lambda I$는 비가역적이다.  [7.4.7](https://render.githubusercontent.com/view/ipynb?commit=a3e483536003d0454458fd57da8d665d19aeca34&enc_url=68747470733a2f2f7261772e67697468756275736572636f6e74656e742e636f6d2f457863656c73696f72434a482f436f64696e675468654d61747269782f613365343833353336303033643034353434353866643537646138643636356431396165636133342f4368617030372532302d25323044696d656e73696f6e2f4368617030372d44696d656e73696f6e2e6970796e62&nwo=ExcelsiorCJH%2FCodingTheMatrix&path=Chap07+-+Dimension%2FChap07-Dimension.ipynb&repository_id=125392345&repository_type=Repository#7.4.7-%ED%96%89%EB%A0%AC%EC%9D%98-%EA%B0%80%EC%97%AD%EC%84%B1)에 의하면 $(A - \lambda I)^T$ 또한 비가역적이다. $(A - \lambda I)^T = A^T - \lambda I$ 이므로, $\lambda$는 $A^T$의 고유값이다.



### 13.3.1 유사성과 대각화 가능성 - Diagonalizability

***Definition*** : 가역행렬 S에 대해 $S^{-1}AS = B$ 가 만족되면 두 정방행렬 $A$와 $B$는 '유사' 또는 '닮은(similar)' 행렬이라고 한다.

***Proposition*** : 유사행렬(similar matrix)들은 동일한 고유값을 가진다.

- **Proof** : $\lambda$ 를 행렬 $A$의 고유값이라 하고, $v$ 를 고유벡터라고 하면, $Av = \lambda v$ 가 성립한다. $S^{-1}AS = B$라 하고, $w = S^{-1}v$라 하면 다음이 성립한다.

$$
\begin{eqnarray} Bw & = & { S }^{ -1 }ASw \\  & = & { S }^{ -1 }AS{ S }^{ -1 }v \\  & = & { S }^{ -1 }Av \\  & = & { S }^{ -1 } \lambda v \\  & = & \lambda S^{-1} v  \\  & = & \lambda w \end{eqnarray}
$$

- 따라서, $\lambda$ 는 행렬 $B$의 고유값이다.

***Example 13.3.11*** : 뒤에서 다루겠지만, 행렬 $U$는 상삼각행렬로, 행렬 $A=\begin{bmatrix} 6 & 3 & -9 \\ 0 & 9 & 15 \\ 0 & 0 & 15 \end{bmatrix}$ 의 고유값은 대각원소들인 $6, 9, 15$ 이다. 행렬 $B=\begin{bmatrix} 92 & -32 & -15 \\ -64 & 34 & 39 \\ 176 & -68 & -99 \end{bmatrix}$ 는 $B=S^{-1}AS$ 인 성질을 가진다. 여기서, $S=\begin{bmatrix} -2 & 1 & 4 \\ 1 & -2 & 1 \\ -4 & 3 & 5 \end{bmatrix}$ 이다. 따라서, $B$의 고유값 또한 $6, 9, 15$ 이다.



***Definition*** : 만약 어떤 정방행렬 $A$가 대각행렬과 유사행렬이면, 즉 대각행렬 $\Lambda$에 대해 $S^{-1}AS = \Lambda$를 만족하는 가역행렬 $S$가 있으면, $A$는 *'대각화가 가능하다(diagonalizable)'* 라고 한다. 