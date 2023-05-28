CHANGES IN VERSION 1.12.1 (2023-05-28):
--------------------------------------
 + `weights_var2test_condi` is enforced to `FALSE` 
 for the permutation test
 

CHANGES IN VERSION 1.11.1 (2023-03-16):
--------------------------------------
 + add parallel support on Windows


CHANGES IN VERSION 1.9.4 (2022-07-20):
--------------------------------------
 + use `scattermore::` in `plot_weights()`


CHANGES IN VERSION 1.8.3 (2022-07-13):
--------------------------------------
 + fixed a bug in NA handling with `which_weights` != "none"


CHANGES IN VERSION 1.8.1 (2022-04-28):
--------------------------------------
 + fixed a small bug in permutation p-values introduced with the 1.8.0 release


CHANGES IN VERSION 1.7.2 (2022-04-13):
--------------------------------------
 + reintroduced `CompQuadForm::davies()` as a longer alternative when the
 `"saddlepoint"` method in `survey::pchisqsum()` fails for computing quadratic 
 form asymptotic p-values
 
 
CHANGES IN VERSION 1.7.1 (2022-02-07):
--------------------------------------
 + added a spaghetti plot functionality


CHANGES IN VERSION 1.5.1 (2021-09-01):
--------------------------------------
 + switching from `CompQuadForm::davies()` to the `"saddlepoint"` method from 
 `survey::pchisqsum()` for computing quadratic form asymptotic p-values


CHANGES IN VERSION 1.3.0 (2021-03-09):
--------------------------------------
 + adding plotting capabilities


CHANGES IN VERSION 1.2.0 (2020-05-20):
--------------------------------------
 + adding adaptive permutations


INITIAL RELEASE - VERSION 1.0.0 (2019-11-26):
---------------------------------------------
 + this is a reboot of the `tcgsaseq` package

