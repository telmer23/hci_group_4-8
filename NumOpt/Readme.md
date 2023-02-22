# Numerical Optimisation
## 28. Give the definition of a descent direction. Draw a simple example explain- ing the properties of a descent direction

In optimization, a descent direction is a direction in which the objective function decreases. More formally, a direction $\mathbf{d}$ is a descent direction at a point $\mathbf{x}$ if there exists a positive scalar $\alpha$ such that $f(\mathbf{x} + \alpha \mathbf{d}) < f(\mathbf{x})$ for some small enough $\alpha > 0$.

To illustrate the properties of a descent direction, consider the following two-dimensional example:

$$f(x,y) = x^2 + y^2$$

The contour plot of this function is a set of concentric circles centered at the origin. Any direction that points towards the origin is a descent direction, as moving in that direction decreases the value of $f$. For example, the vector $\mathbf{d} = (-1,-1)$ is a descent direction at the point $(1,1)$, since $f(1,1) = 2$ and $f(0,0) = 0$.

On the other hand, any direction that points away from the origin is not a descent direction, as moving in that direction increases the value of $f$. For example, the vector $\mathbf{d} = (1,1)$ is not a descent direction at the point $(1,1)$, since $f(1,1) = 2$ and $f(2,2) = 8$.

It is worth noting that not all descent directions are equally good. In particular, the steepest descent direction is the direction of greatest decrease of the objective function, i.e., the negative of the gradient. In the above example, the steepest descent direction at $(1,1)$ is $\nabla f(1,1) = (2,2)$, which points towards the origin. The steepest descent direction is important for many optimization algorithms, as it can lead to faster convergence towards the minimum.

---
## 29. Give the general form of a descent method and show that $d^k = −D^k \nabla f(x^k)$ with $D^k$ symmetric and positive definite is a descent direction. Give three different standard choices for descent directions based on choosing the scaling matrix Dk and discuss their numerical performance.

A descent method is an iterative optimization algorithm that generates a sequence of points ${\mathbf{x}_k}$ that converges to a local minimum of a given objective function $f(\mathbf{x})$. At each iteration $k$, a descent direction $\mathbf{d}k$ is computed, which is a direction that decreases the value of the objective function. The new iterate $\mathbf{x}{k+1}$ is then obtained by taking a step in the descent direction:

$$\mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{d}_k$$

where $\alpha_k$ is a step size parameter that controls the size of the step.

To show that $ \mathbf{d}_k = -D_k \nabla f(\mathbf{x}_k)$ is a descent direction, we need to show that there exists a positive scalar $\alpha_k$ such that $f(\mathbf{x}_k + \alpha_k \mathbf{d}_k) < f(\mathbf{x}_k)$ for some small enough $\alpha_k > 0$. Using a first-order Taylor approximation, we have:

$$f(\mathbf{x}_k + \alpha_k \mathbf{d}_k) \approx f(\mathbf{x}_k) + \alpha_k \nabla f(\mathbf{x}_k)^T \mathbf{d}_k$$

Substituting $\mathbf{d}_k$ and simplifying, we get:

$$f(\mathbf{x}_k + \alpha_k \mathbf{d}_k) \approx f(\mathbf{x}_k) - \alpha_k \nabla f(\mathbf{x}_k)^T D_k \nabla f(\mathbf{x}_k)$$

Since $D_k$ is symmetric and positive definite, its eigenvalues are all positive. Therefore, we can choose $\alpha_k$ small enough such that $\alpha_k \nabla f(\mathbf{x}_k)^T D_k \nabla f(\mathbf{x}_k) > 0$, which implies that $f(\mathbf{x}_k + \alpha_k \mathbf{d}_k) < f(\mathbf{x}_k)$ and thus $\mathbf{d}_k$ is a descent direction.

There are several standard choices for the scaling matrix $D_k$, which can affect the convergence rate and numerical stability of the optimization algorithm:

Steepest Descent Method: $D_k = \mathbf{I}$, the identity matrix. This corresponds to the steepest descent direction, which is the direction of greatest decrease of the objective function. While this method is easy to implement, it can be slow to converge for highly anisotropic problems.
Newton's Method: $D_k = \nabla^2 f(\mathbf{x}_k)$, the Hessian matrix. This corresponds to the exact Newton direction, which takes into account second-order information about the objective function. While this method can converge much faster than the steepest descent method, it can be computationally expensive to compute and may be numerically unstable if the Hessian is poorly conditioned.
Quasi-Newton Methods: $D_k$ is an approximation to the Hessian matrix that is updated iteratively based on the gradient information. The most common example is the Broyden-Fletcher-Goldfarb-Shanno (BFGS) method, which uses a rank-two update formula to update an estimate of the inverse Hessian matrix. This method can have faster convergence than the steepest descent method and can be more stable

---
## 30. Explain three methods how to select the step size in gradient methods.

In gradient methods, selecting an appropriate step size is important for convergence to a local minimum. There are several methods to choose the step size:

Fixed step size: In this method, the step size $\alpha_k$ is kept constant for all iterations. While this method is simple to implement, it can converge slowly for functions with highly anisotropic contours or if the step size is too small or too large.
Backtracking line search: In this method, the step size is computed by starting with a large initial guess for the step size, and then iteratively decreasing it until a sufficient decrease condition is met. The sufficient decrease condition requires that the new objective function value is smaller than the old value by a fixed proportion times the product of the gradient and the step size. Backtracking line search is computationally more expensive than fixed step size, but it can lead to faster convergence.
Exact line search: In this method, the step size is computed by solving a one-dimensional optimization problem that minimizes the objective function along a search direction. This method guarantees that the step size is optimal, but it can be computationally expensive to solve the one-dimensional optimization problem at each iteration.
The choice of step size method depends on the specific optimization problem and the properties of the objective function. Fixed step size can be a good choice for simple problems or when the objective function is well-behaved. Backtracking line search can be a good choice for more complicated problems or when the step size needs to be adapted to the local geometry of the objective function. Exact line search can be a good choice when the objective function is highly non-linear and the other methods are not effective.

---
## 31. Prove the formula for exact line search based on a descent direction d for quadratic functions of the form 
### $$min \frac{1}{2} x^T Ax+b^T x+c$$
To find the exact line search step size, we want to minimize the quadratic function $f(x_k + \alpha d_k)$ along the direction $d_k$, where $\alpha$ is the step size. Using the quadratic function given in the question, we have:

$$\begin{aligned} f(x_k + \alpha d_k) &= \frac{1}{2}(x_k + \alpha d_k)^T A (x_k + \alpha d_k) + b^T (x_k + \alpha d_k) + c \ &= \frac{1}{2}x_k^T A x_k + \alpha d_k^T A x_k + \frac{1}{2} \alpha^2 d_k^T A d_k + b^T x_k + \alpha b^T d_k + c \end{aligned}$$

Taking the derivative of $f(x_k + \alpha d_k)$ with respect to $\alpha$, we get:

$$\frac{\partial f}{\partial \alpha} = d_k^T A x_k + \alpha d_k^T A d_k + b^T d_k$$

Setting this derivative to zero and solving for $\alpha$, we obtain:

$$\alpha = -\frac{d_k^T A x_k + b^T d_k}{d_k^T A d_k}$$

which is the formula for the exact line search step size for quadratic functions with a descent direction $d_k$.

---
## 32. Explain the Armijo step size rule and draw a figure. Prove the sufficient decrease condition and show its relation to the Armijo step size condition.

The Armijo step size rule is a common method for determining the step size in gradient descent methods. It is based on a sufficient decrease condition that ensures that the step size is not too large, so that the objective function decreases sufficiently.

The Armijo rule works by first choosing an initial step size $\alpha$, and then iteratively decreasing it until a sufficient decrease condition is met. The sufficient decrease condition requires that the new objective function value is smaller than the old value by a fixed proportion times the product of the gradient and the step size. Specifically, the Armijo condition is:

$$f(x_k + \alpha \Delta x_k) \leq f(x_k) + c_1 \alpha \nabla f(x_k)^T \Delta x_k$$

where $c_1 \in (0, 1)$ is a fixed constant that controls the amount of decrease required.

The Armijo step size rule can be visualized as follows:

Armijo rule figure

In this figure, the black line represents the objective function $f(x)$. The blue line represents the tangent line at $x_k$, which is used to approximate the objective function locally. The red line represents the next iterate $x_{k+1} = x_k + \alpha \Delta x_k$, which is obtained by taking a step in the direction of the negative gradient.

To prove the sufficient decrease condition, we first note that the Taylor series expansion of the objective function around $x_k$ is:

$$f(x_k + \alpha \Delta x_k) \approx f(x_k) + \alpha \nabla f(x_k)^T \Delta x_k$$

This approximation assumes that the objective function is sufficiently smooth around $x_k$, so that the first-order Taylor approximation is accurate.

Substituting this approximation into the Armijo condition yields:

$$f(x_k) + \alpha \nabla f(x_k)^T \Delta x_k \leq f(x_k) + c_1 \alpha \nabla f(x_k)^T \Delta x_k$$

Simplifying this inequality and rearranging terms gives:

$$f(x_k + \alpha \Delta x_k) \leq f(x_k) + c_1 \alpha \nabla f(x_k)^T \Delta x_k$$

which is the sufficient decrease condition.

The Armijo step size rule ensures that the step size is not too large by iteratively decreasing it until the Armijo condition is satisfied. The larger the value of $c_1$, the more strict the Armijo condition becomes, and the smaller the step size that will be chosen. The Armijo condition can be seen as a way of balancing the need to take large steps (to converge quickly) with the need to take small steps (to ensure that the objective function decreases sufficiently).

---
## 33. What is the “zig-zag” effect? Prove that the differences of successive iterates are orthogonal to each other when using an exact line search.

The "zig-zag" effect refers to the phenomenon where the iterates generated by a gradient descent method oscillate back and forth across the direction of the steepest descent. This can occur when the step size is too large, causing the method to overshoot the minimum and then oscillate back and forth as it tries to converge.

When using an exact line search (i.e., a line search that finds the step size that minimizes the objective function along the search direction), it is possible to show that the differences of successive iterates are orthogonal to each other. This property is known as conjugate gradient descent.

To prove this, we first define the search direction $d_k$ as the negative gradient at iteration $k$, i.e., $d_k = -\nabla f(x_k)$. We then define the step size $\alpha_k$ as the value that minimizes the objective function along the search direction, i.e., $\alpha_k = \arg\min_{\alpha>0} f(x_k + \alpha d_k)$.

The next iterate $x_{k+1}$ is then given by:

$$x_{k+1} = x_k + \alpha_k d_k$$

We can compute the difference between successive iterates as:

$$x_{k+1} - x_k = \alpha_k d_k$$

We want to show that the differences of successive iterates are orthogonal to each other, i.e., that $(x_{k+1} - x_k)^T (x_{k+2} - x_{k+1}) = 0$.

Substituting the expressions for $x_{k+1}$ and $x_k$ yields:

$$(\alpha_k d_k)^T (\alpha_{k+1} d_{k+1}) = 0$$

We know that $d_{k+1} = -\nabla f(x_{k+1})$, so we can rewrite this as:

$$(\alpha_k d_k)^T (-\alpha_{k+1} \nabla f(x_{k+1})) = 0$$

We can simplify this expression by using the definition of the step size:

$$(\alpha_k d_k)^T (-\nabla f(x_k + \alpha_k d_k)) = 0$$

We know that $d_k$ is the negative gradient at iteration $k$, so we can substitute this into the expression:

$$(\alpha_k d_k)^T (-\nabla f(x_k) - \alpha_k \nabla^2 f(x_k) d_k) = 0$$

We can simplify this expression by using the fact that $\nabla f(x_k) = -d_k$:

$$(\alpha_k d_k)^T (d_k - \alpha_k \nabla^2 f(x_k) d_k) = 0$$

We can factor out $\alpha_k$:

$$\alpha_k (d_k^T d_k - \alpha_k d_k^T \nabla^2 f(x_k) d_k) = 0$$

We know that $d_k$ is the steepest descent direction at iteration $k$, so it is orthogonal to all previous search directions $d_0, \ldots, d_{k-1}$. This means that $d_k^T \nabla^2 f(x_j) d_i = 0$ for $j < k$ and $i \leq j$, since the Hessian is symmetric.

Using this fact, we can simplify the expression further:

$$\alpha_k (d_k^T d_k - 0) = 0$$

This shows that the differences of successive iterates are orthogonal to

---
## 34. What is the rate of convergence of the gradient method for quadratic functions when using exact line search. What is the condition number?

The rate of convergence of the gradient method for quadratic functions when using exact line search is linear.

Let $f(x) = \frac{1}{2}x^TQx - b^Tx$ be a quadratic function where $Q$ is a symmetric positive definite matrix and $b$ is a vector. Then the gradient of $f$ at $x_k$ is $\nabla f(x_k) = Qx_k - b$, and the exact line search chooses the step size $\alpha_k$ that minimizes $f(x_k - \alpha_k\nabla f(x_k))$.

The update formula for the gradient method with exact line search is therefore:

$$x_{k+1} = x_k - \alpha_k(Qx_k - b)$$

<!-- Using the fact that $Q$ is symmetric and positive definite, we can write its eigenvalue decomposition as $Q = V\Lambda V^{-1}$, where $V$ is an orthonormal matrix and $\Lambda$ is a diagonal matrix of eigenvalues. Let $y_k = V^{-1}(x_k - x^)$, where $x^$ is the minimizer of $f$. Then we have: -->
Using the fact that $Q$ is symmetric and positive definite, we can write its eigenvalue decomposition as $Q = V\Lambda V^{-1}$, where $V$ is an orthonormal matrix and $\Lambda$ is a diagonal matrix of eigenvalues. Let $y_k = V^{-1}(x_k - x)$, where $x$ is the minimizer of $f$. Then we have:

$$y_{k+1} = y_k - \alpha_kV^{-1}QV y_k = (I - \alpha_k\Lambda)y_k$$

Since $I - \alpha_k\Lambda$ is a diagonal matrix with entries $(1 - \alpha_k\lambda_i)$, where $\lambda_i$ are the eigenvalues of $Q$, the eigenvalues of the matrix $I - \alpha_k\Lambda$ are $1 - \alpha_k\lambda_i$. Therefore, the $i$th component of $y_{k+1}$ is $(1 - \alpha_k\lambda_i)y_{k,i}$, and we have:

$$|y_{k+1}|^2 = \sum_{i=1}^n (1 - \alpha_k\lambda_i)^2 y_{k,i}^2 \leq \lambda_{\max}^2 \sum_{i=1}^n (1 - \alpha_k\lambda_i)^2 y_{k,i}^2$$

where $\lambda_{\max}$ is the maximum eigenvalue of $Q$.

Taking the minimum of the summands inside the sum gives:

$$\sum_{i=1}^n (1 - \alpha_k\lambda_i)^2 y_{k,i}^2 \geq \frac{1}{n}\left(\sum_{i=1}^n (1 - \alpha_k\lambda_i)y_{k,i}\right)^2 = \frac{1}{n}|y_k - y_{k+1}|^2$$

Using the fact that $y_k$ converges to the minimizer $x^*$, we can choose $\alpha_k$ such that $\alpha_k\lambda_{\max} < 2$, which gives:

$$|y_{k+1}|^2 \leq (1 - \frac{2\alpha_k\lambda_{\max}}{n})|y_k|^2$$

This shows that the rate of convergence of the gradient method with exact line search for quadratic functions is linear, with a convergence factor of $1 - \frac{2\alpha_k\lambda_{\max}}{n}$.

The condition number of $Q$ is $\kappa = \frac{\lambda_{\max}}{\lambda_{\min}}$, where $\lambda_{\min}$ is the minimum eigenvalue of $Q$. A large condition number indicates that $Q$ is ill-conditioned, meaning

---
## 35. What is a Lipschitz continuous gradient and show how it is related to the norm of the Hessian matrix? What is the Lipschitz constant of the gradient of a linear least squares problem of the form


---
## 36. Prove the descent lemma for a differentiable function f with Lipschitz- continuous gradient. Interpret the inequality of the descent lemma in terms of an upper bound to the function f. Use the descent lemma to also prove the sufficient decrease lemma.
Descent Lemma:
Let $f:\mathbb{R}^n\to\mathbb{R}$ be a differentiable function with Lipschitz-continuous gradient and let $\mathbf{x}, \mathbf{y} \in \mathbb{R}^n$ such that $\mathbf{y} = \mathbf{x} - t\nabla f(\mathbf{x})$ for some $t > 0$. Then, we have the following inequality:
$$f(\mathbf{y}) \leq f(\mathbf{x}) - \frac{t}{2} |\nabla f(\mathbf{x})|_2^2.$$

Proof:
Since $f$ is differentiable, we can apply the mean value theorem to obtain

$$ $$

for some $\mathbf{z}$ on the line segment connecting $\mathbf{x}$ and $\mathbf{y}$. Using the Lipschitz continuity of $\nabla f$, we have

$$ $$

where $L$ is the Lipschitz constant of $\nabla f$. Since $\nabla f(\mathbf{y}) = \nabla f(\mathbf{x}) - t\nabla f(\mathbf{x})$, we have

$$ $$

Substituting these expressions into the mean value theorem, we obtain

$$ f(y) = f(\mathbf{x}) - t \nabla f(x)^T \nabla f(\mathbf{x}) - t \nabla f(\mathbf{x})^T(\nabla f(\mathbf{x}) - t \nabla f(\mathbf{x})) ≤ f(\mathbf{x})- \frac{t}{2}||\nabla f(\mathbf{x})||_2^2   $$

Interpretation of the Descent Lemma:
The inequality of the descent lemma states that the value of the objective function $f$ decreases by at least $\frac{t}{2}|\nabla f(\mathbf{x})|_2^2$ when we move from $\mathbf{x}$ to $\mathbf{y}$ along the direction of the negative gradient. In other words, the descent lemma provides an upper bound on the amount that the objective function can decrease by when moving along a descent direction. This bound is proportional to the magnitude of the gradient of $f$ at $\mathbf{x}$.

Sufficient Decrease Lemma:
Let $f:\mathbb{R}^n\to\mathbb{R}$ be a differentiable function with Lipschitz-continuous gradient and let $

---

## 38. Give the general form of non-linear least-squares problems. Show how the Gauss-Newton method is obtained from performing a first-order Taylor approximation. What is the Levenberg-Marquadt method?

The general form of a non-linear least-squares problem is:

$$min_{x} f(x) = \sum_{i=1}^{m} [y_{i} - f_{i}(x)]^{2}$$

where $x$ is the parameter vector to be estimated, $y_{i}$ is the observed data, $f_{i}(x)$ is the model prediction for the $i$-th observation, and $m$ is the number of observations.

To solve this problem, we can use the Gauss-Newton method, which is obtained by performing a first-order Taylor approximation of the non-linear model around the current estimate of $x$. The Taylor approximation is given by:

$$f_{i}(x + \Delta x) \approx f_{i}(x) + J_{i}\Delta x$$

where $J_{i}$ is the Jacobian matrix of $f_{i}(x)$ evaluated at $x$. Substituting this approximation into the objective function and neglecting higher-order terms, we obtain the linearized least-squares problem:

$$min_{\Delta x} ||J\Delta x + r||^{2}$$

where $J$ is the Jacobian matrix of the model function $f(x)$ evaluated at the current estimate $x$, and $r$ is the residual vector given by $r_{i} = y_{i} - f_{i}(x)$. This linear least-squares problem can be solved using standard methods such as the QR decomposition or the singular value decomposition.

The Gauss-Newton method then updates the estimate of $x$ as:

$$x_{k+1} = x_{k} + \Delta x$$

where $\Delta x$ is the solution of the linear least-squares problem.

However, the Gauss-Newton method may encounter convergence issues when the Jacobian matrix is ill-conditioned or the objective function has a non-convex shape. To address these issues, the Levenberg-Marquardt method introduces a damping factor $\lambda$ to the Gauss-Newton update rule:

$$(J^{T}J + \lambda I)\Delta x = J^{T}r$$

where $I$ is the identity matrix. When $\lambda$ is large, the method behaves like a gradient descent method and moves slowly towards the minimum. When $\lambda$ is small, the method behaves like the Gauss-Newton method and converges faster. The value of $\lambda$ is adjusted at each iteration based on the reduction in the objective function.
