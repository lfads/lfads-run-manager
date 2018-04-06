# LFADS Hyperparameters

The following is an nearly exhaustive list of hyperparameters that affect the training and posterior mean sampling of the model. If the parameter name begins with `c_`, this parameter will be passed to the LFADS Python code directly, without the `c_` prefix. All of these values are specified in the `RunParams` instance that accompanies each `Run`.

While there are many parameters, you will likely care about only a small subset of them. We have color-coded the hyperparameters below according to how often they require tuning:

<style>
th { font-weight: bold; }
tr.hp-common td:first-child { background: #F44336; font-weight: bold; color: #ffffff; }
tr.hp-medium td:first-child { background: #FFCDD2; }
tr.hp-rare td:first-child { background: #ffffff; color: #455A64; }

tbody.hp tr td:first-child { font-family: "Roboto Mono","Courier New",Courier,monospace; }

</style>

<table>
<thead>
<th>Tuning Frequency</th>
<th>Description</th>
</thead>
<tbody>
<tr class="hp-common">
<td>Common</td>
<td>Typically requires tuning on a per-project basis and/or is important to set appropriately upfront.</td>
</tr>
<tr class="hp-medium">
<td>Occasional</td>
<td>Might be adjusted for fine tuning.</td>
</tr>
<tr class="hp-rare">
<td>Rare</td>
<td>Infrequently requires tuning and/r primarily intended for advanced users.</td>
</tr>
</table>

## Run Manager logistics and data Processing
<table>
<thead>
<th>Name</th>
<th>Default</th>
<th>Description</th>
</thead>
<tbody class="hp">

<tr class="hp-medium">
<td>name</td>
<td><code>''</code></td>
<td>Name of this set of parameters, used for convenience only. <em>Note that this parameter does not affect either the param or data hash.</em></td>
</tr>

<tr class="hp-rare">
<td>version</td>
<td>n/a</td>
<td><strong>This value you should not assign directly</strong>, as it will automatically be set to match the version of the <code>RunCollection</code> to which it is added. This is used for graceful backwards compatibility. <em>Note that this parameter does not affect either the param or data hash.</em></td>
</tr>

<tr class="hp-common">
<td>spikeBinMs</td>
<td>2</td>
<td>Spike bin width in milliseconds. This must be an integer multiple of the original bin width provided by the `Run` class by `generateCountsForDataset`.</td>
</tr>

</tbody>
</table>

## TensorFlow logistics
<table>
<thead>
<th>Name</th>
<th>Default</th>
<th>Description</th>
</thead>
<tbody class="hp">

<tr class="hp-medium">
<td>c_allow_gpu_growth</td>
<td>true</td>
<td>Allow the GPU to dynamically allocate memory instead of allocating all the GPU's memory at the start</td>
</tr>

<tr class="hp-rare">
<td>c_max_ckpt_to_keep</td>
<td>5</td>
<td>Max number of checkpoints to keep (rolling)</td>
</tr>

<tr class="hp-rare">
<td>c_max_ckpt_to_keep_lve</td>
<td>5</td>
<td>Max number of checkpoints to keep for lowest validation error models (rolling)</td>
</tr>

<tr class="hp-rare">
<td>c_device</td>
<td><code>'gpu:0'</code></td>
<td>Which visible GPU or CPU to use. Note that GPUs are typically scheduled by setting `CUDA_VISIBLE_DEVICES` rather than using this parameter.</td>
</tr>

</tbody>
</table>

## Optimization

Rather put the learning rate on an exponentially decreasiong schedule,
the current algorithm pays attention to the learning rate, and if it
isn't regularly decreasing, it will decrease the learning rate.  So far,
it works fine, though it is not perfect.

<table>
<thead>
<th>Name</th>
<th>Default</th>
<th>Description</th>
</thead>
<tbody class="hp">

<tr class="hp-rare">
<td>c_learning_rate_init</td>
<td>0.01</td>
<td>Initial learning rate</td>
</tr>

<tr class="hp-medium">
<td>c_learning_rate_decay_factor</td>
<td>0.95</td>
<td>Factor by which to decrease the learning rate if progress isn't being made.</td>
</tr>

<tr class="hp-rare">
<td>c_learning_rate_n_to_compare</td>
<td>6</td>
<td>Number of previous costs current cost has to be worse than, in order to lower learning rate.</td>
</tr>

<tr class="hp-rare">
<td>c_learning_rate_stop</td>
<td>0.00001</td>
<td>Stop training when the learning rate reaches this threshold.</td>
</tr>

<tr class="hp-rare">
<td>c_max_grad_norm</td>
<td>200</td>
<td>Max norm of gradient before clipping. This sets a value, above which, the gradients will be clipped.  This hp
is extremely useful to avoid an infrequent, but highly pathological
problem whereby the gradient is so large that it destroys the
optimization by setting parameters too large, leading to a vicious cycle
that ends in NaNs.  If it's too large, it's useless, if it's too small,
it essentially becomes the learning rate. It's pretty insensitive, though.</td>
</tr>

<tr class="hp-rare">
<td>trainToTestRatio</td>
<td>4</td>
<td>Ratio of training vs testing trials used.</td>
</tr>

<tr class="hp-common">
<td>c_batch_size</td>
<td>256</td>
<td>Number of trials to use during each training pass. The total trial count must be &ge; <code>c_batch_size * (trainToTestRatio + 1)</code>.</td>
</tr>

<tr class="hp-rare">
<td>c_cell_clip_value</td>
<td>5</td>
<td>Max value recurrent cell can take before being clipped. If your optimizations start "NaN-ing out", reduce this value so that the values of the network don't grow out of control.  Typically, once this parameter is set to a reasonable value, one stops having numerical problems.</td>
</tr>

</tbody>
</table>

## Overfitting

If controller is heavily penalized, then it won't have any output. If dynamics are heavily penalized, then generator won't make dynamics.  Note this l2 penalty is only on the recurrent portion of the RNNs, as dropout is also available, penalizing the feed-forward connections.

<table>
<thead>
<th>Name</th>
<th>Default</th>
<th>Description</th>
</thead>
<tbody class="hp">

<tr class="hp-rare">
<td>c_temporal_spike_jitter_width</td>
<td>0</td>
<td>Enables jittering spike times during training. It appears that the system will happily fit spikes (blessing or curse, depending).  You may not want this.  Jittering the spikes a bit may help (-/+ bin size, as specified here). The idea is to prevent LFADS from trying to learn very fine temporal structure in the data if you believe this to be noise.</td>
</tr>

<tr class="hp-medium">
<td>c_keep_prob</td>
<td>0.95</td>
<td>Fraction of units to randomly drop during each training pass. Dropout is done on the input data, on controller inputs (from encoder), and on outputs from generator to factors.</td>
</tr>

<tr class="hp-medium">
<td>c_l2_gen_scale</td>
<td>500</td>
<td>L2 regularization cost for the generator only.</td>
</tr>

<tr class="hp-rare">
<td>c_l2_con_scale</td>
<td>500</td>
<td>L2 regularization cost for the controller only.</td>
</tr>

<tr class="hp-rare">
<td>c_co_mean_corr_scale</td>
<td>0</td>
<td>Cost of correlation (through time) in the means of controller output.</td>
</tr>

</tbody>
</table>

## Underfitting

If the primary task of LFADS is "filtering" of data and not
generation, then it is possible that the KL penalty is too strong.
Empirically, we have found this to be the case.  So we add a
hyperparameter in front of the the two KL terms (one for the initial
conditions to the generator, the other for the controller outputs).
You should always think of the the default values as 1.0, and that
leads to a standard VAE formulation whereby the numbers that are
optimized are a lower-bound on the log-likelihood of the data. When
these 2 HPs deviate from 1.0, one cannot make any statement about
what those LL lower bounds mean anymore, and they cannot be compared
(AFAIK).

Sometimes the task can be sufficiently hard to learn that the
optimizer takes the 'easy route', and simply minimizes the KL
divergence, setting it to near zero, and the optimization gets
stuck. The same possibility is true for the L2 regularizer. One wants a simple generator, for scientific
reasons, but not at the expense of hosing the optimization. The last 5 parameters help avoid that by by getting the
optimization to 'latch' on to the main optimization, and only turning on the regularizers gradually by increasing their weighting in the overall cost functions later.

<table>
<thead>
<th>Name</th>
<th>Default</th>
<th>Description</th>
</thead>
<tbody class="hp">

<tr class="hp-medium">
<td>c_kl_ic_weight</td>
<td>1</td>
<td>Strength of KL weight on initial conditions KL penalty.</td>
</tr>

<tr class="hp-medium">
<td>c_kl_co_weight</td>
<td>1</td>
<td>Strength of KL weight on controller output KL penalty.</td>
</tr>

<tr class="hp-medium">
<td>c_kl_start_step</td>
<td>0</td>
<td>Start increasing KL weight after this many steps.</td>
</tr>

<tr class="hp-medium">
<td>c_kl_increase_steps</td>
<td>900</td>
<td>Number of steps over which the KL weight increases.</td>
</tr>

<tr class="hp-rare">
<td>c_l2_start_step</td>
<td>0</td>
<td>Start increasing L2 weight after this many steps.</td>
</tr>

<tr class="hp-medium">
<td>c_l2_increase_steps</td>
<td>900</td>
<td>Number of steps over which the L2 weight increases.</td>
</tr>

<tr class="hp-rare">
<td>c_l2_start_step</td>
<td>0</td>
<td>Start increasing L2 weight after this many steps</td>
</tr>

<tr class="hp-medium">
<td>scaleIncreaseStepsWithDatasets</td>
<td>true</td>
<td>If true, <code>c_kl_increase_steps</code> and <code>c_l2_increase_steps</code> will be multiplied by the number of datasets in a stitching run.</td>
</tr>

</tbody>
</table>

## External inputs

If there are observed inputs, there are two ways to add that observed
input to the model.  The first is by treating as something to be
inferred, and thus encoding the observed input via the encoders, and then
input to the generator via the "inferred inputs" channel.  Second, one
can input the input directly into the generator.  This has the downside
of making the generation process strictly dependent on knowing the
observed input for any generated trial.

<table>
<thead>
<th>Name</th>
<th>Default</th>
<th>Description</th>
</thead>
<tbody class="hp">

<tr class="hp-common">
<td>c_ext_input_dim</td>
<td>0</td>
<td>Number of external, known (or observed) inputs.</td>
</tr>

<tr class="hp-medium">
<td>c_inject_ext_input_to_gen</td>
<td>false</td>
<td>Should the known inputs be input to model via encoders (false) or injected directly into generator (true)?</td>
</tr>

</tbody>
</table>

## Controller and inferred inputs

The controller will be more powerful if it can see the encoding of the entire
trial.  However, this allows the controller to create inferred inputs that are
acausal with respect to the actual data generation process.  For example, the data
generator could have an input at time `t`, but the controller, after seeing the
entirety of the trial could infer that the input is coming a little before
time `t`, because there are no restrictions on the data the controller sees.
One can force the controller to be causal (with respect to perturbations in
the data generator) so that it only sees forward encodings of the data at time
`t` that originate at times before or at time `t`.  One can also control the data
the controller sees by using an input lag (forward encoding at time `t-tlag`
for controller input at time `t`.  The same can be done in the reverse direction
(controller input at time `t` from reverse encoding at time `t+tlag`, in the
case of an acausal controller).  Setting this lag > 0 (even lag=1) can be a
powerful way of avoiding very spiky decodes. Finally, one can manually control
whether the factors at time `t-1` are fed to the controller at time `t`.

If you don't care about any of this, and just want to smooth your data, set
`do_causal_controller = False`, `do_feed_factors_to_controller = True`, `controller_input_lag = 0`.

<table>
<thead>
<th>Name</th>
<th>Default</th>
<th>Description</th>
</thead>
<tbody class="hp">

<tr class="hp-common">
<td>c_co_dim</td>
<td>4</td>
<td>Number of inferred inputs (controller outputs). This parameter critically controls whether or not there is a controller (along with controller encoders placed into the LFADS graph. If equal to 0, no controller will be added.</td>
</tr>

<tr class="hp-rare">
<td>c_prior_ar_atau</td>
<td>10</td>
<td>Initial autocorrelation of AR(1) priors (in time bins)</td>
</tr>

<tr class="hp-rare">
<td>c_do_train_prior_ar_atau</td>
<td>true</td>
<td>Is the value for atau an initial value (true) or the constant value (false)?</td>
</tr>

<tr class="hp-rare">
<td>c_prior_ar_nvar</td>
<td>0.1</td>
<td>Initial noise variance for AR(1) priors</td>
</tr>

<tr class="hp-rare">
<td>c_do_train_prior_ar_nvar</td>
<td>true</td>
<td>Is the value for the noise var an initial value (true) or the constant value (false)?</td>
</tr>

<tr class="hp-medium">
<td>c_do_causal_controller</td>
<td>false</td>
<td>Restrict input encoder from seeing the future?</td>
</tr>

<tr class="hp-medium">
<td>c_do_feed_factors_to_controller</td>
<td>true</td>
<td>Should <code>factors[t-1]</code> be input to controller at time t? Strictly speaking, feeding either the factors or the rates to the controller violates causality, since the g0 gets to see all the data. This may or may not be only a theoretical concern.</td>
</tr>

<tr class="hp-rare">
<td>c_feedback_factors_or_rates</td>
<td><code>'factors'</code></td>
<td>Feedback the factors or the rates to the controller? Set to either <code>'factors'</code> or <code>'rates'</code></td>
</tr>

<tr class="hp-medium">
<td>c_controller_input_lag</td>
<td>1</td>
<td>Time lag on the encoding to controller <code>t-lag</code> for forward, <code>t+lag</code> for reverse.</td>
</tr>

<tr class="hp-medium">
<td>c_ci_enc_dim</td>
<td>128</td>
<td>Network size for controller input encoder.</td>
</tr>

<tr class="hp-medium">
<td>c_con_dim</td>
<td>128</td>
<td>Controller dimensionality.</td>
</tr>

<tr class="hp-rare">
<td>c_co_prior_var_scale</td>
<td>0.1</td>
<td>Variance of control input prior distribution.</td>
</tr>

</tbody>
</table>

# Encoder and initial conditions for generator

Note that the dimension of the initial conditions is separated from the
dimensions of the generator initial conditions (and a linear matrix will
adapt the shapes if necessary).  This is just another way to control
complexity.  In all likelihood, setting the IC dims to the size of the
generator hidden state is just fine.

For the initial condition prior variance parameters, it's best to leave them alone.
IThe defaults should be fine for most cases, irregardless of other parameters.
If you don't want the prior variance to be learned, set the
following values to the same thing: `ic_prior_var_min`, `ic_prior_var_scale`, `ic_prior_var_max`.
The prior mean will still be learned. If you really want to limit the information from encoder to decoder,
increase `ic_post_var_min` above 0.

<table>
<thead>
<th>Name</th>
<th>Default</th>
<th>Description</th>
</thead>
<tbody class="hp">

<tr class="hp-rare">
<td>c_num_steps_for_gen_ic</td>
<td><code>MAXINT</code></td>
<td>Number of steps to train the generator initial condition.</td>
</tr>

<tr class="hp-common">
<td>c_ic_dim</td>
<td>64</td>
<td>Dimensionality of the initial conditions.</td>
</tr>

<tr class="hp-rare">
<td>c_ic_enc_dim</td>
<td>128</td>
<td>Network size for IC encoder.</td>
</tr>

<tr class="hp-rare">
<td>c_ic_prior_var_min</td>
<td>0.1</td>
<td>Minimum variance of IC prior distribution</td>
</tr>

<tr class="hp-rare">
<td>c_ic_prior_var_scale</td>
<td>0.1</td>
<td>Variance of IC prior distribution</td>
</tr>

<tr class="hp-rare">
<td>c_ic_prior_var_max</td>
<td>0.1</td>
<td>Maximum variance of IC prior distribution</td>
</tr>

<tr class="hp-rare">
<td>c_ic_post_var_min</td>
<td>0.0001</td>
<td>Minimum variance of IC posterior distribution</td>
</tr>

</tbody>
</table>

## Generator network, factors, rates

Controlling the size of the generator is one way to control complexity of
the dynamics (there is also l2, which will squeeze out unnecessary
dynamics also).  The modern deep learning approach is to make these cells
as large as tolerable (from a waiting perspective), and then regularize
them to death with drop out or whatever. It is not clear if this is correct
for the LFADS application or not.

<table>
<thead>
<th>Name</th>
<th>Default</th>
<th>Description</th>
</thead>
<tbody class="hp">

<tr class="hp-rare">
<td>c_cell_weight_scale</td>
<td>1.0</td>
<td>Input scaling for input weights in generator. The combined recurrent and input weights of the encoder and controller cells are by default set to scale at ws/sqrt(#inputs) with ws=1.0.  You can change this scaling with this parameter.</td>
</tr>

<tr class="hp-common">
<td>c_gen_dim</td>
<td>100</td>
<td>Generator network size/td>
</tr>

<tr class="hp-rare">
<td>c_gen_cell_input_weight_scale</td>
<td>1.0</td>
<td>Input scaling for input weights in generator, which will be divided by sqrt(#inputs)</td>
</tr>

<tr class="hp-rare">
<td>c_gen_cell_rec_weight_scale</td>
<td>1.0</td>
<td>Input scaling for recurrent weights in generator.</td>
</tr>

<tr class="hp-common">
<td>c_factors_dim</td>
<td>50</td>
<td>Dimensionality of factors read out from generator network. This provides dimensionality reduction from generator dimensionality down to factors and then back out to the neural rates. <em>Note that this property does affect the data and param hashes, unlikely the other <code>c_</code> prefixed parameters, which only affect the param hash.</td>
</tr>

<tr class="hp-rare">
<td>c_output_dist</td>
<td>'poisson'</td>
<td>Type of output distribution for rates, either <code>'poisson'</code> or <code>'gaussian'</code></td>
</tr>

</tbody>
</table>

## Stitching multi-session models

<table>
<thead>
<th>Name</th>
<th>Default</th>
<th>Description</th>
</thead>
<tbody class="hp">

<tr class="hp-common">
<td>c_do_train_readin</td>
<td>true</td>
<td>For stitching models, make the readin matrices trainable (true) or fix them to equal the alignment matrices (false). The per-session readin matrices map from neurons to input factors which are fed into the shared encoder. These are initialized by the alignment matrices and can subsequently be fixed or made trainable.</td>
</tr>

<tr class="hp-common">
<td>useAlignmentMatrix</td>
<td>false</td>
<td>Whether to use an alignment matrix when stitching datasets together./td>
</tr>

<tr class="hp-rare">
<td>useSingleDatasetAlignmentMatrix</td>
<td>false</td>
<td>When only using a single dataset, it is also possible to use a readin matrix that reduces the dimensionality of the spikes before inputting these input factors to the encoder networks. If set true, this will set up this readin matrix and seed it with an alignment matrix computed using PCA.</td>
</tr>

</tbody>
</table>

## Posterior sampling

<table>
<thead>
<th>Name</th>
<th>Default</th>
<th>Description</th>
</thead>
<tbody class="hp">

<tr class="hp-common">
<td>posterior_mean_kind</td>
<td><code>'posterior_sample_and_average'</code></td>
<td>Mechanism to obtain the posterior mean. Either <code>'posterior_sample_and_average'</code> to take a specified number of samples from the posterior distribution, run them through the model, and average the results. Or <code>'posterior_push_mean'</code> to use the posterior mean of the ICs and inputs and push those through the model directly. Since there are nonlinearities in the network, this need not be equivalent to the mean of the samples, but in practice it's usually pretty close, and is much faster to compute. <em>Note that this parameter does not affect either the param or data hash.</em></td>
</tr>

<tr class="hp-medium">
<td>num_samples_posterior</td>
<td>512</td>
<td>Number of samples of the posterior to use when using <code>'posterior_sample_and_average'</code>. <em>Note that this parameter does not affect either the param or data hash.</em></td>
</tr>

</tbody>
</table>
