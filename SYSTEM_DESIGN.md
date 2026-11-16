# Linear and nonlinear system design

Current model specification: 9 November 2026.

We use a synthetic two-variable system to study how chronology uncertainty affects causal inference. X evolves on its own. Y retains part of its previous state, receives a delayed contribution from X, and experiences its own random disturbances. The linear and nonlinear cases share the same memory times, background variability and coupling gain. They differ in the function that transforms X into the input to Y.

The main parameter choices are a 20-year time step, a 500-year intrinsic memory for X, a 200-year intrinsic memory for Y, and a long-run coupling gain of 4. These choices determine the discrete coefficients. In particular, the coefficients 0.904837 and 0.380650 are calculated from the specified memory and gain; they are not two separately chosen tuning parameters.

This specification describes the process generator and the curves in [the model-truth figure](figures/02_check_causal_pipeline/00_model_truth.png). It follows [the generator](toolbox/process.py), the settings and generator calls in [the main experiment](02_check_causal_pipeline.py), and the current saved [baseline parameters](data/processed/02_check_causal_pipeline/baseline_parameters.json). The saved parameter file is a local generated artifact; the source files contain the definitions needed to reproduce the model.

## 1. Time, units and notation

The subscript t denotes a discrete step advancing forward in physical time. It is not a BP-age coordinate. All time parameters in the Python code use kyr, so the current step and imposed delay are

$$
\Delta t=0.020\ \mathrm{kyr}=20\ \mathrm{yr},
\qquad d=0.020\ \mathrm{kyr}=20\ \mathrm{yr}.
$$

The number of delayed steps is

$$
\ell=\frac{d}{\Delta t}=1.
$$

The generator requires this ratio to be a positive integer. Thus, X at the preceding step enters the equation for Y at the current step. The symbol ell here denotes the process delay; it is unrelated to a CCM library length.

Within the retained record, elapsed time increases as age decreases:

$$
\mathrm{elapsed\ time}_t=t\Delta t,
\qquad
\mathrm{age}_t=30\ \mathrm{kyr}-t\Delta t,
$$

$$
\mathrm{calendar\ year}_t=1950-1000\,\mathrm{age}_t,
$$

where age is expressed in kyr BP. Plotting older ages on the left does not reverse the causal direction.

| Symbol | Meaning | Current value or role |
|---|---|---|
| X, Y | Synthetic process values | Arbitrary units |
| T_X, T_Y | Intrinsic decay times | 500 and 200 years |
| a_X, a_Y | Fractions of the previous state retained per step | Derived from memory and step length |
| S_X | Theoretical stationary SD of X | 3.202563 |
| S_Y | Stationary SD of Y when coupling is absent | 1.666667 |
| u_t, v_t | New process disturbances at step t | Independent, zero-mean Gaussian innovations |
| sigma_u, sigma_v | SDs of those innovations | Derived from background SDs and memory |
| z_t | X expressed in units of S_X | X_t / S_X |
| g(X) | Normalized driver transformation | Linear or saturating nonlinear function |
| c | Coupling coefficient applied at each step | 0.380650 for both coupled cases |
| K | Long-run mean response to a held unit g(X) | 4 for the coupled cases; 0 for the independent case |

Because g is dimensionless, c and K carry the output scale of Y. Their numerical values depend on the chosen arbitrary units; they are not universal dimensionless measures of causal strength.

## 2. The complete process equations

Both coupled systems use

$$
X_t=a_X X_{t-1}+u_t,
$$

$$
Y_t=a_Y Y_{t-1}+c\,g(X_{t-\ell})+v_t.
$$

The innovation distributions are

$$
u_t\overset{\mathrm{iid}}{\sim}\mathcal N(0,\sigma_u^2),
\qquad
v_t\overset{\mathrm{iid}}{\sim}\mathcal N(0,\sigma_v^2).
$$

The complete u and v sequences are independent of one another. Each innovation is also independent of earlier process states. The resulting X and Y sequences are temporally correlated because their equations retain previous states.

There is no Y term in the X equation. The coupled systems therefore have the prescribed cross-variable direction X to Y only. Their self-memory terms remain part of the dynamics.

The independent control sets c=0, giving

$$
X_t=a_X X_{t-1}+u_t,
\qquad
Y_t=a_Y Y_{t-1}+v_t.
$$

Process noise remains present in every case. The main experiment's zero measurement noise means that no additional measurement disturbance is added when observing these already stochastic processes.

## 3. Deriving the coefficients from intrinsic memory

An intrinsic memory time describes how fast a perturbation fades when the subsequent forcing and innovations are held the same. Consider two Y trajectories with identical future inputs but different initial Y values. Their difference obeys

$$
\delta Y_t=a_Y\delta Y_{t-1},
\qquad
\delta Y_{t+m}=a_Y^m\delta Y_t.
$$

We specify an exponential decay over physical time:

$$
\delta Y_{t+m}
=\exp\left(-\frac{m\Delta t}{T_Y}\right)\delta Y_t.
$$

Matching the two expressions gives

$$
\boxed{a_Y=\exp\left(-\frac{\Delta t}{T_Y}\right)},
\qquad
\boxed{a_X=\exp\left(-\frac{\Delta t}{T_X}\right)}.
$$

The inverse conversion, for the positive coefficients used here, is

$$
T_Y=-\frac{\Delta t}{\ln a_Y},
\qquad
T_X=-\frac{\Delta t}{\ln a_X}.
$$

For the current choices,

$$
a_X=e^{-20/500}\approx0.9607894392,
\qquad
a_Y=e^{-20/200}\approx0.9048374180.
$$

Y retains about 90.48% of an existing perturbation after one 20-year step. After 200 years, it retains exp(-1), or about 36.79%. The 200-year memory is an e-folding time, not a cutoff after which all earlier influence disappears.

The half-life is a different quantity:

$$
t_{1/2}=T\ln 2.
$$

It is approximately 346.57 years for X and 138.63 years for Y. The positive coefficients below one give stable, monotonic decay of an initial perturbation.

For an unforced stationary AR(1) process, the autocorrelation at a separation of m steps has the same exponential form:

$$
\rho_X(m)=a_X^{|m|}
=\exp\left(-\frac{|m|\Delta t}{T_X}\right),
$$

$$
\rho_{Y,\,c=0}(m)=a_Y^{|m|}
=\exp\left(-\frac{|m|\Delta t}{T_Y}\right).
$$

The second expression applies to Y without coupling. The total autocorrelation of driven Y also reflects the history of X, as explained below.

## 4. Deriving innovation amplitudes while preserving background variability

Changing memory without adjusting the innovations would also change the stationary process variance. We avoid that confounding by specifying background SDs first and deriving the innovation SDs.

For stationary X, the innovation is independent of the preceding state, so

$$
\operatorname{Var}(X_t)
=a_X^2\operatorname{Var}(X_{t-1})+\sigma_u^2.
$$

Setting both process variances to S_X squared gives

$$
S_X^2=a_X^2S_X^2+\sigma_u^2,
$$

and hence

$$
\boxed{\sigma_u=S_X\sqrt{1-a_X^2}}.
$$

The same calculation for Y with c=0 gives

$$
\boxed{\sigma_v=S_Y\sqrt{1-a_Y^2}}.
$$

The currently specified background scales are

$$
S_X=\frac{1}{\sqrt{1-0.95^2}}\approx3.2025630761,
$$

$$
S_Y=\frac{1}{\sqrt{1-0.8^2}}\approx1.6666666667.
$$

The constants 0.95 and 0.8 in these amplitude definitions are not the current a_X and a_Y. The current AR coefficients are derived from the 500-year and 200-year memories at the 20-year step.

Substituting the current coefficients gives

$$
\sigma_u\approx0.8880035436,
\qquad
\sigma_v\approx0.7095954382.
$$

S_Y is the background SD without coupling. It is not a constraint that forces the total SD of driven Y to remain 1.666667. Coupling adds variability with its own temporal structure.

## 5. The linear and nonlinear driver functions

The standardized driver is

$$
z_t=\frac{X_t}{S_X}.
$$

X has zero theoretical mean and a stationary Gaussian distribution. Therefore its standardized marginal distribution is standard normal:

$$
Z\sim\mathcal N(0,1).
$$

Successive z values remain temporally correlated. Standardizing the amplitude does not remove memory.

### Linear coupling

The linear case uses

$$
\boxed{g_{\mathrm{lin}}(X)=z=\frac{X}{S_X}}.
$$

Thus the driver has theoretical mean zero and SD one:

$$
\mathbb E[g_{\mathrm{lin}}(X)]=0,
\qquad
\operatorname{Var}[g_{\mathrm{lin}}(X)]=1.
$$

A doubling of X doubles the direct input to Y. The per-step coefficient multiplying raw X is c/S_X, approximately 0.118858, rather than c itself, because the model uses standardized X.

### Nonlinear coupling

The nonlinear case uses a normalized hyperbolic tangent:

$$
\boxed{g_{\mathrm{nonlin}}(X)=\frac{\tanh(\beta z)}{s_\beta}},
\qquad \beta=2.
$$

The hyperbolic tangent is

$$
\tanh(q)=\frac{e^q-e^{-q}}{e^q+e^{-q}}.
$$

It passes smoothly through zero, increases monotonically, and approaches -1 or +1 for sufficiently negative or positive input. The positive constant beta controls how quickly the transition occurs. The current beta=2 is a prescribed shape parameter.

The normalization is defined using the theoretical standard normal distribution:

$$
s_\beta=\sqrt{\operatorname{Var}[\tanh(\beta Z)]}.
$$

Because tanh is odd and the normal distribution is symmetric,

$$
\mathbb E[\tanh(\beta Z)]=0,
$$

so

$$
s_\beta^2
=\frac{1}{\sqrt{2\pi}}
\int_{-\infty}^{\infty}
\tanh^2(\beta z)e^{-z^2/2}\,dz.
$$

The code evaluates this fixed integral with 256-point Gauss-Hermite quadrature. If x_i and w_i are the quadrature nodes and weights for the weight exp(-x squared), the implemented approximation is

$$
s_\beta^2\approx\frac{1}{\sqrt{\pi}}
\sum_{i=1}^{256}w_i\tanh^2(\beta\sqrt{2}\,x_i).
$$

For beta=2, the implementation uses

$$
s_2\approx0.7970327690,
\qquad
\boxed{g_{\mathrm{nonlin}}(X)
=\frac{\tanh(2X/S_X)}{0.7970327690}}.
$$

The theoretical normalization gives zero mean and unit variance, up to numerical quadrature accuracy. This constant is calculated from the theoretical distribution and is fixed across realizations. The code does not estimate it separately from each simulated record.

Normalization makes the marginal variance of the direct forcing c g(X) equal to c squared in both coupled cases, up to the quadrature approximation:

$$
\operatorname{Var}[c\,g(X)]\approx c^2.
$$

For the linear transformation the equality is exact. Matching this marginal scale does not make the two driver functions have identical temporal covariance, and it does not force the total Y variances to match.

The nonlinear function's slope with respect to z is

$$
\frac{d g_{\mathrm{nonlin}}}{dz}
=\frac{\beta}{s_\beta}
\left[1-\tanh^2(\beta z)\right].
$$

It is largest at zero and approaches zero as the response saturates. Positive X gives positive forcing, and negative X gives negative forcing. The model has no abrupt threshold, hysteresis or reverse causal link.

## 6. Deriving c from the long-run gain K

We define K through a simple experiment: hold X at a fixed value x, keep its transformed value g(x) constant, and let the mean response of Y settle. Write that limiting mean as mu_Y(x). Because v has zero mean,

$$
\mu_Y(x)=a_Y\mu_Y(x)+c\,g(x).
$$

Solving gives

$$
\mu_Y(x)=\frac{c}{1-a_Y}g(x).
$$

We therefore define

$$
\boxed{K=\frac{c}{1-a_Y}},
\qquad
\boxed{c=(1-a_Y)K}.
$$

For the current K=4,

$$
c=\left(1-e^{-20/200}\right)\times4\approx0.3806503279.
$$

The Y equation can now be written as

$$
Y_t=a_Y Y_{t-1}+(1-a_Y)K\,g(X_{t-\ell})+v_t,
$$

or, equivalently,

$$
\boxed{
Y_t-Y_{t-1}
=(1-a_Y)\left[K\,g(X_{t-\ell})-Y_{t-1}\right]+v_t.
}
$$

This form is easy to interpret. At each step, the deterministic update closes about 9.52% of the gap between the preceding Y value and the current target K g(X). A random innovation is added as well. The target changes as X changes.

If the target is held constant after the imposed delay has elapsed, the mean after m updates is

$$
\mathbb E[Y_m]
=K g(x)+a_Y^m\left(\mathbb E[Y_0]-K g(x)\right).
$$

Starting from mean zero, this becomes

$$
\mathbb E[Y_m]=K g(x)(1-a_Y^m)
=K g(x)\left[1-e^{-m\Delta t/T_Y}\right].
$$

After 200 years of this constant input, the mean has completed about 63.21% of the adjustment; after 600 years, about 95.02%. The response begins after the 20-year equation delay and then continues to adjust over the 200-year intrinsic timescale. These are different roles for time.

## 7. The current equations with numerical values

The following rounded equations show the actual parameterization. Calculations use the unrounded values recorded in the saved metadata.

The common X process is

$$
X_t=0.9607894392\,X_{t-1}+u_t,
\qquad
u_t\sim\mathcal N\left(0,(0.8880035436)^2\right).
$$

The linear Y process is

$$
\boxed{
Y_t=0.9048374180\,Y_{t-1}
+0.3806503279\frac{X_{t-1}}{3.2025630761}+v_t.
}
$$

The nonlinear Y process is

$$
\boxed{
Y_t=0.9048374180\,Y_{t-1}
+0.3806503279
\frac{\tanh\left(2X_{t-1}/3.2025630761\right)}{0.7970327690}
+v_t.
}
$$

For both,

$$
v_t\sim\mathcal N\left(0,(0.7095954382)^2\right).
$$

The independent control retains the same memory coefficients and innovation amplitudes but removes the entire X-driven term from the Y equation.

## 8. Reading the curves in the model-truth figure

The lower panels show the limiting mean of Y after holding X fixed. They are not plots of the noisy generated values, and they are not fitted regressions of observed Y on observed X.

Their horizontal coordinate is z=x/S_X. Their vertical coordinate is K g(x), although the panel titles display g(x). For K=4, the two curves are

$$
\boxed{\mu_{Y,\mathrm{lin}}(z)=4z},
$$

$$
\boxed{\mu_{Y,\mathrm{nonlin}}(z)
=\frac{4\tanh(2z)}{0.7970327690}}.
$$

| Held z | Linear mean response | Nonlinear mean response |
|---:|---:|---:|
| 0 | 0.00 | 0.00 |
| 0.5 | 2.00 | 3.82 |
| 1 | 4.00 | 4.84 |
| 2 | 8.00 | 5.02 |
| 3 | 12.00 | 5.02 |

The responses at negative z are equal in magnitude and opposite in sign. The nonlinear mean approaches finite plateaus:

$$
\lim_{z\to\pm\infty}\mu_{Y,\mathrm{nonlin}}(z)
=\pm\frac{K}{s_2}
\approx\pm5.018614.
$$

Its slope with respect to the plotted z coordinate is

$$
\frac{d\mu_{Y,\mathrm{nonlin}}}{dz}
=\frac{2K}{s_2}\left[1-\tanh^2(2z)\right].
$$

At z=0 this is approximately 10.0372, whereas the linear slope is 4 everywhere. Thus K=4 does not mean that every local slope of the nonlinear curve equals 4. For slopes with respect to raw X, divide these slopes by S_X.

These plateaus bound the mean response to a held driver. They do not bound individual noisy Y values, because Gaussian process disturbances can produce excursions beyond them. In the independent control the held-driver mean is zero, while Y still fluctuates because v remains present.

## 9. Why the intrinsic memory of Y differs from its total persistence

Repeated substitution in the Y equation gives the following stationary-limit representation, after the initial transient has decayed:

$$
Y_t
=c\sum_{j=0}^{\infty}a_Y^j g(X_{t-\ell-j})
+\sum_{j=0}^{\infty}a_Y^j v_{t-j}.
$$

The first term is the accumulated contribution of earlier X values. The second is the accumulated contribution of Y's own disturbances. Using c=(1-a_Y)K, the first term can be written as

$$
Y_t^{\mathrm{forced}}
=K\sum_{j=0}^{\infty}w_j g(X_{t-\ell-j}),
\qquad
w_j=(1-a_Y)a_Y^j,
$$

with

$$
\sum_{j=0}^{\infty}w_j=1.
$$

The forced component of Y is therefore a weighted average of past transformed X, multiplied by K. The weights decay on Y's intrinsic timescale, but the X values being averaged are themselves persistent. X has a 500-year memory, and its nonlinear transformation can further change the covariance of the forcing.

This explains why estimating an autocorrelation time from the full driven Y record need not return 200 years. The 200-year value describes Y's own relaxation under matched subsequent inputs. It does not summarize every source of persistence in the coupled record.

Similarly, the prescribed 20-year direct delay does not require a cross-correlation peak or a CCM lag-scan optimum to occur at exactly 20 years. Such summaries depend on the process histories and on the estimator as well as on the direct delay.

## 10. Which values are chosen, and which are derived?

| Model choice | Current specification | Resulting calculation |
|---|---|---|
| Native step | Delta t = 20 yr | Determines the duration represented by one update |
| X intrinsic memory | T_X = 500 yr | a_X = exp(-Delta t / T_X) |
| Y intrinsic memory | T_Y = 200 yr | a_Y = exp(-Delta t / T_Y) |
| Background amplitudes | S_X = 3.202563; S_Y = 1.666667 | Innovation SDs = background SD times sqrt(1 - a squared) |
| Direct delay | d = 20 yr | ell = d / Delta t = 1 |
| Long-run gain | K = 4 | c = (1 - a_Y) K |
| Nonlinear shape | beta = 2 | tanh(beta z), with its theoretical normalization |
| Warm-up multiplier | Five times the longest memory | 2.5 kyr, plus the separate required lag history |

The algebra explains the derived coefficients, but it does not establish that 500 years, 200 years, K=4 or beta=2 are uniquely correct physical values. They define the present synthetic benchmark. The current specification does not provide an empirical calibration of these values to a particular climate system.

The memory choices make X evolve more slowly than Y's intrinsic relaxation. The gain sets the mean response to a sustained standardized forcing. For example, a held g(x)=1 gives a mean response of 4, compared with an unforced Y background SD of about 1.67. This is a description of the benchmark amplitude, not a calculation of the signal-to-noise ratio of the full coupled record.

The existing [coupling-strength experiment](02_1_check_coupling_strength.py) compares K=1, 2 and 4 using three paired process seeds while fixing the other process settings. It examines how native method scores respond to gain. It does not establish that K=4 is optimal, and it does not test the choice of a 200-year memory. Any use of that comparison must also retain the estimator settings recorded for its saved results.

K is a parameter of the generating model. CCM correlations, GPDC statistics, information-flow scores and Granger scores have their own definitions; none is required to equal K.

## 11. Changing memory or step length without changing the wrong quantity

If T_Y changes while the intended long-run gain and unforced background SD remain fixed, three quantities must be recalculated together:

$$
a_Y'=e^{-\Delta t/T_Y'},
\qquad
c'=(1-a_Y')K,
\qquad
\sigma_v'=S_Y\sqrt{1-(a_Y')^2}.
$$

At the current 20-year step, K=4 and S_Y=1.666667 give:

| T_Y (yr) | a_Y | c | Innovation SD sigma_v |
|---:|---:|---:|---:|
| 100 | 0.818731 | 0.725077 | 0.956963 |
| 200 | 0.904837 | 0.380650 | 0.709595 |
| 400 | 0.951229 | 0.195082 | 0.514141 |

A longer memory retains more of each preceding state. It therefore needs a smaller per-step forcing coefficient to reach the same constant-input mean, and smaller innovations to preserve the same unforced background SD. These rows are algebraic illustrations, not additional completed simulations.

If c were instead held fixed while a_Y changed, the long-run gain would change to

$$
K'=\frac{c}{1-a_Y'}.
$$

If the innovation SD were held fixed, the unforced stationary SD would change to

$$
S_Y'=\frac{\sigma_v}{\sqrt{1-(a_Y')^2}}.
$$

Both effects would mix a memory change with another change in the experiment. Preserving K and S_Y still does not preserve the total driven Y variance, because the filtering of X changes with memory.

Changing the time step requires the same recalculation for both variables:

$$
a_X'=e^{-\Delta t'/T_X},
\qquad
a_Y'=e^{-\Delta t'/T_Y},
$$

$$
\sigma_u'=S_X\sqrt{1-(a_X')^2},
\qquad
\sigma_v'=S_Y\sqrt{1-(a_Y')^2},
\qquad
c'=(1-a_Y')K.
$$

The imposed delay must also remain representable:

$$
\ell'=\frac{d}{\Delta t'}\in\{1,2,3,\ldots\}.
$$

The generator receives discrete coefficients explicitly. Changing its dt argument alone does not recalculate them; the calling analysis script performs these conversions.

These conversions preserve the specified intrinsic exponential decay, unforced variance and constant-input gain. They do not make simulations at different steps identical subsamples of one common trajectory. In particular, the delayed coupled update is a specified discrete model; the project does not claim it is an exact discretization of a joint continuous-time stochastic differential equation.

## 12. Initialization, warm-up and matched comparisons

X is initialized from its stationary marginal distribution. Here, index zero denotes the start of generation before warm-up; the elapsed-time index in Section 1 is reset at the first retained sample. In the code,

$$
X_0=\frac{u_0}{\sqrt{1-a_X^2}},
$$

which gives X_0 distributed as N(0, S_X squared). Y is initialized at zero and is allowed to warm up under the evolving driver and its innovations.

The warm-up length is

$$
N_{\mathrm{burn}}
=\left\lceil\frac{5\max(T_X,T_Y)}{\Delta t}\right\rceil
=125,
\qquad
N_{\mathrm{burn}}\Delta t=2.5\ \mathrm{kyr}.
$$

The implementation allocates and discards the required ell-step lag history separately, so the first retained array index is

$$
i_{\mathrm{start}}=N_{\mathrm{burn}}+\ell.
$$

A finite warm-up reduces the influence of initialization; it does not guarantee that Y has reached its exact stationary distribution. At the current T_Y, an initial Y difference subjected to identical later inputs decays over the 2.5-kyr warm-up by

$$
\exp(-2500/200)=\exp(-12.5)\approx3.73\times10^{-6}.
$$

The retained process contains both endpoints of a 30-kyr span:

$$
N_{\mathrm{retained}}=\frac{30}{0.020}+1=1501.
$$

The main causal analysis uses the common 29-to-1 kyr BP window:

$$
N_{\mathrm{analysis}}=\frac{29-1}{0.020}+1=1401.
$$

For a given main-experiment replicate, the three process forms reuse the same process seed and innovations. Their X trajectory is therefore identical, and their Y disturbances are paired; the driver transformation or absence of coupling changes the Y response. This comparison holds incidental random variation fixed. The three forms are not three independent realizations of one model.

Irregular observations and chronology ensembles are constructed after generating these processes. They affect the observations or assigned times, rather than changing T_X, T_Y, K, g or the equation delay. Ten chronology members therefore do not represent ten independently generated X/Y histories.

## 13. Connection to the implementation

| Mathematical quantity or operation | Code location or name |
|---|---|
| Delta t, T_X, T_Y, K, d | `DT`, `MEMORY_X`, `MEMORY_Y`, `COUPLING_GAIN`, `TRUE_LAG` in the root analysis script |
| S_X, S_Y | `BACKGROUND_X`, `BACKGROUND_Y` |
| Derivation of AR coefficients and innovation SDs | Generator preparation in `02_check_causal_pipeline.py` |
| Per-step coupling c | `coupling=(1 - ar_y) * gain` passed to `generate_process` |
| Standardized driver z | `driver = x / driver_std` in `toolbox/process.py` |
| Nonlinear normalization s_2 | `coupling_transform_std`, calculated with 256 Hermite nodes |
| Nonlinear transformation | `np.tanh(2 * driver) / coupling_transform_std` |
| Dynamic update of Y | The loop using `ar_y`, delayed `driver`, `coupling` and `noise_y` |
| Held-driver mean K g(X) | `plot_process_truth` in `02_check_causal_pipeline.py` |

The default arguments in the general `generate_process` function are not the settings of the main experiment. The main script passes its derived coefficients explicitly. When describing a saved experiment, use its recorded parameters and the corresponding generator calls rather than assuming the function defaults were used.
