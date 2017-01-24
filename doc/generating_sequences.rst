.. _generating_sequences:

Generating Random Sequences
###########################

Throughout the code, anytime the polymer needs to be labelled with a binary
sequence, a markov chain (whose states are the values of the label at each bead)
is used to initialize the chain.

If we call the labels $L_1 = A$ and $L_2 = B$ for now, this Markov chain is
defined in terms of two parameters, the fraction of A and the non-unit
eigenvalue of the transition state matrix,
$$M = \left(\begin{bmatrix} p_{1,1} & p{2,1} \\ p_{1,2} & p{2,2} \end{bmatrix}\right)$$
which is simply $1 - p_{1,2} - p_{2,1}$, where
$p_{i,j} \coloneq P(X^{(n+1)} == L_j | X^{(n)} == L_i)$. Of course, $p_{1,2} = 1
- p_{1,1}$ and $p_{2,1} = 1 - p_{2,2}$, so there are just two parameters of the
transition matrix, $p_1 = p_{1,1}$ and $p_2 = p_{2,2}$.

Straightforward calculation gives that the average run length of A's is
$\frac{1}{1 - p_1}$, that the fraction of A's is just $\frac{1 - p_2}{2 - p_1 -
p_2}$, and that the eigenvalues of the transition matrix are 1 and $-1 + p_1 +
p_2$.

So to specify a particular Markov chain defined by $p_1$ and $p_2$, use the
quantites $\texttt{lam} = -1 + p_1 + p_2$ and $\texttt{fA} = \frac{1 - p_2}{2 - p_1 -
p_2}$ as simulation inputs.
