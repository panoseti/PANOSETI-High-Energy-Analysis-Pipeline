# Introduction to robust estimation

Here we introduce the concept of robust estimation, which is a statistical method that aims to provide reliable parameter estimates even in the presence of outliers or deviations from model assumptions. We motivate the basic approach of robust fitting with M-estimators starting from the fitting of a Gaussian model to a set of measurements using maximum likelihood estimation (MLE). A more generic approach would not start from the assumption of a Gaussian model, but we will use it here to illustrate the basic ideas, see [Huber (1981)](https://onlinelibrary.wiley.com/doi/book/10.1002/0471725250) for a more general treatment.

### Maximum likelihood estimation of a Gaussian model

Consider the problem of fitting a Gaussian model to a set of measurements ${x_i}$ ($1\le i\le N$) which may contain a significant set of outliers from an different distribution that we do not want to try to model explicitly. We assume the measurements are independent and identically distributed (i.i.d.) and that the Gaussian model has mean $m$ and standard deviation $s$. The probability density function of a single measurement $x_i$ is given by:

$$p(x_i;m,s)=\frac{1}{\sqrt{2\pi}s}e^{-\frac{(x_i-m)^2}{2s^2}}$$

A classic approach to estimating the parameters of the Gaussian is with (unbinned) maximum likelihood. We define the likelihood function for the a set of measurements is the product of these probabilities expressed as a function of $m$ and $s$. 

$$\mathcal L(m,s)=\prod_{i=1}^N \frac{1}{\sqrt{2\pi}s}e^{-\frac{(x_i-m)^2}{2s^2}}$$

We maximize the likelihood function, or alternatively minimize the negative log-likelihood (score) function, which is given by:

$$\ell(m,s)=-\log\mathcal L(m,s)=\sum_{i=1}^N \frac{(x_i-m)^2}{2s^2} + N\log s\tag{1}$$

where we have dropped the constant term $\frac{N}{2}\log(2\pi)$ which does not depend on the parameters. For this model the solution can be derived analytically by taking the derivatives of the score function with respect to $m$ and $s$ and setting them to zero. The derivative with respect to $m$ is given by:

$$\frac{\partial\ell}{\partial m} = \sum_{i=1}^N\frac{x_i-m}{s^2}=0 \iff m = \frac{1}{N}\sum_{i=1}^N x_i$$

and the derivative with respect to $s$ is given by:

$$\frac{\partial\ell}{\partial s} = -\sum_{i=1}^N\frac{(x_i-m)^2}{s^3}+\frac{N}{s}=0 \iff s^2 = \frac{1}{N}\sum_{i=1}^N (x_i-m)^2$$

A maximum likelihood fit of a Gaussian model to a set of data therefore results in estimates of the model parameters $m$ and $s$ given by the sample mean and sample standard deviation, respectively.

The problem should be clear: the sample mean and standard deviation are not robust to outliers. A single outlier can result in an unbounded bias the estimates of $m$ and $s$. 

### Robust estimation and M-estimators

Robust estimation techniques aim to limit the influence of outliers on the parameter estimates. One approach to robust estimation is to use M-estimators, which generalize the maximum likelihood estimation by minimizing a different objective function that reduces the impact of outliers.

To motivate this we can rewrite the score function in equation (1) identically as:

$$\ell(m,s) = \sum_{i=1}^N \rho\left(\frac{x_i-m}{s}\right) + N\log s\tag{2}$$

where $\rho(u)=u^2/2$, is called the loss function. We can rewrite the derivative of the score function with respect to $m$ and $s$ using the function $\psi(u) = \rho'(u)$, which is called the influence function, as:

$$\frac{\partial\ell}{\partial m} = \sum_{i=1}^N \psi\left(\frac{x_i-m}{s}\right)\frac{1}{s}$$

and 

$$\frac{\partial\ell}{\partial s} = -\sum_{i=1}^N \psi\left(\frac{x_i-m}{s}\right)\frac{x_i-m}{s^2} + \frac{N}{s}$$

Setting these derivatives to zero gives two equations that can be solved to find the estimates of $m$ and $s$:

$$\frac{1}{N} \sum_{i=1}^N \psi\left(\frac{x_i-m}{s}\right) = 0\tag{3}$$

and 

$$\frac{1}{N} \sum_{i=1}^N \psi\left(\frac{x_i-m}{s}\right)\frac{x_i-m}{s} = \beta\tag{4}$$

where we define $\beta = 1$. 

This expression is identical to the Gaussian MLE case, and we would get the same estimates of $m$ and $s$ as before. However, we can now choose a different loss function $\rho(u)$ that is less sensitive to outliers. For example, we can use the *Huber loss function*, which is defined as:

$$\rho(u) = \begin{cases} \frac{1}{2}u^2 & |u| \leq \delta \\ \delta|u| - \frac{1}{2}\delta^2 & |u| > \delta \end{cases}$$

where $\delta$ is a tuning parameter that controls the threshold for outlier suppression. This function increases quadratically for small values of $u$ (like the Gaussian case) but linearly for large values of $u$, which reduces the influence of outliers on the parameter estimates.

The corresponding influence function is given by:

$$\psi(u) = \begin{cases} u & |u| \leq \delta \\ \delta \cdot \text{sign}(u) & |u| > \delta \end{cases}$$

which is bounded by $\pm\delta$, limiting the effect that any single outlier can have on the estimate of the parameters. With this choice of loss function, we can solve equations (3) and (4) iteratively or using a minimizer to obtain robust estimates of $m$ and $s$ that are less sensitive to outliers in the data. Within the context of robust estimation, the parameter $m$ is conventionally referred to as the *location* parameter, and $s$ is referred to as the *scale* parameter. The value of $\beta$ must be adjusted to calibrate the scale estimator to be consistent with the Gaussian case.

An alternative way to view equations (3) and (4) is using the concept of *weights*. We can define a weight function:

$$w(u) = \begin{cases} \psi(u)/u\tag{5} & u\ne0 \\ \psi'(0) & u=0 \end{cases}$$

For the Huber loss function we have:

$$w(u) = \begin{cases} 1 & |u| \leq \delta \\ \delta/|u| & |u| > \delta \end{cases}$$

If we write $w_i= w\left(\frac{x_i-m}{s}\right)$, we can then rewrite equations (3) and (4) as:

$$\frac{1}{N} \sum_{i=1}^N w_i\frac{x_i-m}{s} = 0$$

and 

$$\frac{1}{N} \sum_{i=1}^N w_i\left(\frac{x_i-m}{s}\right)^2 = \beta$$

Which result in the following equations for $m$ and $s$:

$$m = \frac{\sum_{i=1}^N w_i x_i}{\sum_{i=1}^N w_i}\tag{6}$$

and

$$s^2 = \frac{\sum_{i=1}^N w_i(x_i-m)^2}{\beta\sum_{i=1}^N w_i}\tag{7}$$

These have the form of weighted mean and weighted variance, where the weights $w_i$ are determined by the influence function (though note that the weights are not independent of $m$ and $s$). The weights effectively downweight the influence of outliers, leading to more robust parameter estimates. These expressions can be used iteractively to find robust estimates of the location and scale parameters. 

