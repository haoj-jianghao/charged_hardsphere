SUBROUTINE menu

	USE menu_parameters

	PRINT ("(1X,/15X, 'Integral Equation'//&

				&15X, '1. Particle-particle OZ'/&
				&15X, '2. Site-site OZ for long chain'/&
				&15X, '3. Site-site OZ for molecules with one type of site'/)")

	WRITE(*,"(15X, A)", ADVANCE = "NO") 'Your choice : '
	
	READ *, model
	
	PRINT ("(1X)")

	PRINT ("(1X,/15X, 'Potential model'//&

				&15X, '1. Lennard-Jones'/&
				&15X, '2. Square-Well'/&
				&15X, '3. Hard-Sphere'/)")

	WRITE(*,"(15X, A)", ADVANCE = "NO") 'Your choice : '
	
	READ *, potential
	
	PRINT ("(1X)")

	PRINT ("(1X,/15X, 'Closure'//&

				&15X, '1. Percus-Yevick [PY]'/&
				&15X, '2. Hypernetted Chain [HNC]'/&
				&15X, '3. Mean Spherical Approximation [MSA]'/&
				&15X, '4. Thermodynamically Consistent [TC]'/)")

	WRITE(*,"(15X, A)", ADVANCE = "NO") 'Your choice : '
	
	READ *, closure
	
	PRINT ("(1X)")

	PRINT ("(1X,/15X, 'Iteration method'//&

				&15X, '1. Variational Method (for RISM)'/&
				&15X, '2. Newton-Raphson+Direct Substitution [Gillan]'/&
				&15X, '3. Variational Method+Direct Substitution with Relaxation'/&
				&15X, '4. Variational Method (for particle-particle OZ)'/)")

	WRITE(*,"(15X, A)", ADVANCE = "NO") 'Your choice : '
	
	READ *, iteration

	PRINT ("(1X)")

END SUBROUTINE menu
