## GSL features to be implemented

- [x] Polynomial Evaluation 
	- [x] real coefficients and variable [`gsl_poly_eval()`]
	- [x] real coefficients, complex variable [`gsl_poly_complex_eval()`]
	- [x] complex coefficients and variable [`gsl_complex_poly_complex_eval()`].
	- [x] Polynomial derivatives evaluation [`gsl_poly_eval_derivs()`].
- [ ] Divided Differences
	- [ ] representation calculation [`gsl_poly_dd_init()`]
	- [ ] evaluation ['gsl_poly_dd_eval()']
	- [ ] conversion to Taylor expansion [`gsl_poly_dd_taylor()`]
	- [ ] Hermite representation calculation [`gsl_poly_dd_hermite_init()`]
	- [ ] Hermite representation evaluation [also `gsl_poly_dd_eval()`]
- [ ] Quadratic Equations
	- [x] Calculation of real roots of quadratic equation [`gsl_poly_solve_quadratic()`]
	- [ ] Calculation of complex roots of quadratic equation [`gsl_poly_complex_solve_quadratic()`]
- [ ] Cubix Equations
	- [ ] Calculation of real roots of cubic equation [`gsl_poly_solve_cubic()`]
	- [ ] Calculation of complex roots of cubic equation [`gsl_poly_complex_solve_cubic()`]
- [ ] General Polynomial Equations
