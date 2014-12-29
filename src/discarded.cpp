



#if 0

      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (csystem=="cubic"){
	if (isym==0){
	  /* Result: E - E0 = 0.5 * V0 * C11 * f^2 */
	  eps1 = f;
	}
	else if (isym==1){
	  /* Result: E - E0 = V0 * C12 * f^2*/
	  eps1 = f;
	  eps2 = f;
	}
	else if (isym==2){
	  /* Result: E - E0 = 2 * V0 * C44 * f^2*/
	  eps4 = f;
	}
	else quit=true;
      }
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      else if (csystem=="hexagonal"){
	if (isym==0){
	  /* Result: E - E0 = 0.5 * V0 * C11 * f^2 */
	  eps1 = f;
	}
	else if (isym==1){
	  /* Result: E - E0 = 0.5 * V0 * C33 * f^2 */
	  eps3 = f;
	}
	else if (isym==2){
	  /* Result: E - E0 = V0 * C12 * f^2 */
	  eps1 = f;
	  eps2 = f;
	}
	else if (isym==3){
	  /* Result: E - E0 = V0 * C13 * f^2 */
	  eps1 = f;
	  eps3 = f;
	}
	else if (isym==4){
	  /* Result: E - E0 = 2 * V0 * C44 * f^2*/
	  eps4 = f;
	}
	else quit=true;
      }
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      else if (csystem=="orthorombic"){
	if (isym==0){
	  /* Result: E - E0 = 0.5 * V0 * C11 * f^2 */
	  eps1 = f;
	}
	else if (isym==1){
	  /* Result: E - E0 = 0.5 * V0 * C22 * f^2 */
	  eps2 = f;
	}
	else if (isym==2){
	  /* Result: E - E0 = 0.5 * V0 * C33 * f^2 */
	  eps3 = f;
	}
	else if (isym==3){
	  /* Result: E - E0 = 2.0 * V0 * C44 * f^2 */
	  eps4 = f;
	}
	else if (isym==4){
	  /* Result: E - E0 = 2.0 * V0 * C55 * f^2 */
	  eps5 = f;
	}
	else if (isym==5){
	  /* Result: E - E0 = 2.0 * V0 * C66 * f^2 */
	  eps6 = f;
	}
	else if (isym==6){
	  /* Result: E - E0 = V0 * C12 * f^2 */
	  eps1 = f;
	  eps2 = f;
	}
	else if (isym==7){
	  /* Result: E - E0 = V0 * C13 * f^2 */
	  eps1 = f;
	  eps3 = f;
	}
	else if (isym==8){
	  /* Result: E - E0 = V0 * C23 * f^2 */
	  eps2 = f;
	  eps3 = f;
	}
	else quit=true;
      }
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      else if (csystem=="triclinic"){
	if (isym==0){
	  /* Result: E - E0 = 0.5 * V0 * C11 * f^2 */
	  eps1 = f;
	}
	else if (isym==1){
	  /* Result: E - E0 = 0.5 * V0 * C22 * f^2 */
	  eps2 = f;
	}
	else if (isym==2){
	  /* Result: E - E0 = 0.5 * V0 * C33 * f^2 */
	  eps3 = f;
	}
	else if (isym==3){
	  /* Result: E - E0 = 2.0 * V0 * C44 * f^2 */
	  eps4 = f;
	}
	else if (isym==4){
	  /* Result: E - E0 = 2.0 * V0 * C55 * f^2 */
	  eps5 = f;
	}
	else if (isym==5){
	  /* Result: E - E0 = 2.0 * V0 * C66 * f^2 */
	  eps6 = f;
	}
	else if (isym==6){
	  /* Result: E - E0 = V0 * C12 * f^2 */
	  eps1 = f;
	  eps2 = f;
	}
	else if (isym==7){
	  /* Result: E - E0 = V0 * C13 * f^2 */
	  eps1 = f;
	  eps3 = f;
	}
	else if (isym==8){
	  /* Result: E - E0 = 2.0 * V0 * C14 * f^2 */
	  eps1 = f;
	  eps4 = f;
	}
	else if (isym==9){
	  /* Result: E - E0 = 2.0 * V0 * C15 * f^2 */
	  eps1 = f;
	  eps5 = f;
	}
	else if (isym==10){
	  /* Result: E - E0 = 2.0 * V0 * C16 * f^2 */
	  eps1 = f;
	  eps6 = f;
	}

	else if (isym==11){
	  /* Result: E - E0 = V0 * C23 * f^2 */
	  eps2 = f;
	  eps3 = f;
	}
	else if (isym==12){
	  /* Result: E - E0 = 2.0 * V0 * C24 * f^2 */
	  eps2 = f;
	  eps4 = f;
	}
	else if (isym==13){
	  /* Result: E - E0 = 2.0 * V0 * C25 * f^2 */
	  eps2 = f;
	  eps5 = f;
	}
	else if (isym==14){
	  /* Result: E - E0 = 2.0 * V0 * C26 * f^2 */
	  eps2 = f;
	  eps6 = f;
	}








	else if (isym==8){
	  /* Result: E - E0 = V0 * C23 * f^2 */
	  eps2 = f;
	  eps3 = f;
	}
	else quit=true;








	if (isym==0){
	  /* Result: E - E0 = 0.5 * V0 * C11 * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==1){
	  /* Result: E - E0 = 0.5 * V0 * C22 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==2){
	  /* Result: E - E0 = 0.5 * V0 * C33 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0 + f;
	}
	else if (isym==3){
	  /* Result: E - E0 = 2 * V0 * C44 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(1,2) = f;
	}
	else if (isym==4){
	  /* Result: E - E0 = 2 * V0 * C55 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(0,2) = f;
	}
	else if (isym==5){
	  /* Result: E - E0 = 2 * V0 * C66 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(0,1) = f;
	}

	else if (isym==6){
	  /* Result: E - E0 = 0.5 * V0 * (C11 + C22 + 2.0 * C12) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==7){
	  /* Result: E - E0 = 0.5 * V0 * (C11 + C33 + 2.0 * C13) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0 + f;
	}
	else if (isym==8){
	  /* Result: E - E0 = 0.5 * V0 * (C11 + 4.0 * C44 + 2.0 * C14) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(1,2) = f;
	}
	else if (isym==9){
	  /* Result: E - E0 = 0.5 * V0 * (C11 + 4.0 * C55 + 2.0 * C15) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(0,2) = f;
	}
	else if (isym==10){
	  /* Result: E - E0 = 0.5 * V0 * (C11 + 4.0 * C66 + 2.0 * C16) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(0,1) = f;
	}


	else if (isym==11){
	  /* Result: E - E0 = 0.5 * V0 * (C22 + C33 + 2.0 * C23) * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0 + f;
	}
	else if (isym==12){
	  /* Result: E - E0 = 0.5 * V0 * (C22 + 4.0 * C44 + 2.0 * C24) * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(1,2) = f;
	}
	else if (isym==13){
	  /* Result: E - E0 = 0.5 * V0 * (C22 + 4.0 * C55 + 2.0 * C25) * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(0,2) = f;
	}
	else if (isym==14){
	  /* Result: E - E0 = 0.5 * V0 * (C22 + 4.0 * C66 + 2.0 * C26) * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(0,1) = f;
	}


	else if (isym==15){
	  /* Result: E - E0 = 0.5 * V0 * (C33 + 4.0 * C44 + 2.0 * C34) * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0 + f;
	  alpha.elem(1,2) = f;
	}
	else if (isym==16){
	  /* Result: E - E0 = 0.5 * V0 * (C33 + 4.0 * C55 + 2.0 * C35) * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0 + f;
	  alpha.elem(0,2) = f;
	}
	else if (isym==17){
	  /* Result: E - E0 = 0.5 * V0 * (C33 + 4.0 * C66 + 2.0 * C36) * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0 + f;
	  alpha.elem(0,1) = f;
	}


	else if (isym==18){
	  /* Result: E - E0 = 0.5 * V0 * (4.0 * C44 + 4.0 * C55 + 8.0 * C45) * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(1,2) = f;
	  alpha.elem(0,2) = f;
	}
	else if (isym==19){
	  /* Result: E - E0 = 0.5 * V0 * (4.0 * C44 + 4.0 * C66 + 8.0 * C46) * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(1,2) = f;
	  alpha.elem(0,1) = f;
	}

	else if (isym==20){
	  /* Result: E - E0 = 0.5 * V0 * (4.0 * C55 + 4.0 * C66 + 8.0 * C56) * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(0,2) = f;
	  alpha.elem(0,1) = f;
	}
	else quit=true;
      }




      else if (csystem=="monoclinic"){
	if (isym==0){
	  /* Result: E - E0 = 0.5 * V0 * C11 * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==1){
	  /* Result: E - E0 = 0.5 * V0 * C22 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==2){
	  /* Result: E - E0 = 0.5 * V0 * C33 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0 + f;
	}
	else if (isym==3){
	  /* Result: E - E0 = 2 * V0 * C44 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(1,2) = f;
	}
	else if (isym==4){
	  /* Result: E - E0 = 2 * V0 * C55 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(0,2) = f;
	}
	else if (isym==5){
	  /* Result: E - E0 = 2 * V0 * C66 * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0;
	  alpha.elem(0,1) = f;
	}

	else if (isym==6){
	  /* Result: E - E0 = 0.5 * V0 * (C11 + C22 + 2.0 * C12) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0;
	}
	else if (isym==7){
	  /* Result: E - E0 = 0.5 * V0 * (C11 + C33 + 2.0 * C13) * f^2 */
	  alpha.elem(0,0) = 1.0 + f;
	  alpha.elem(1,1) = 1.0;
	  alpha.elem(2,2) = 1.0 + f;
	}

	else if (isym==8){
	  /* Result: E - E0 = 0.5 * V0 * (C22 + C33 + 2.0 * C23) * f^2 */
	  alpha.elem(0,0) = 1.0;
	  alpha.elem(1,1) = 1.0 + f;
	  alpha.elem(2,2) = 1.0 + f;
	}
	// Monoclinic 1: C15, C25, C35, C46 neq. 0
	// Monoclinic 2: C16, C26, C36, C45 neq. 0
	else if (isym==9){
	  if (csystem_sub==0){
	    /* Result: E - E0 = 0.5 * V0 * (C11 + 4.0 * C55 + 2.0 * C15) * f^2 */
	    alpha.elem(0,0) = 1.0 + f;
	    alpha.elem(1,1) = 1.0;
	    alpha.elem(2,2) = 1.0;
	    alpha.elem(0,2) = f;
	  }
	  else {

	  }
	}
	else if (isym==10){
	  if (csystem_sub==0){	  
	    /* Result: E - E0 = 0.5 * V0 * (C22 + 4.0 * C55 + 2.0 * C25) * f^2 */
	    alpha.elem(0,0) = 1.0;
	    alpha.elem(1,1) = 1.0 + f;
	    alpha.elem(2,2) = 1.0;
	    alpha.elem(0,2) = f;
	  }
	  else {

	  }
	}
	else if (isym==11){
	  if (csystem_sub==0){
	    /* Result: E - E0 = 0.5 * V0 * (C33 + 4.0 * C55 + 2.0 * C35) * f^2 */
	    alpha.elem(0,0) = 1.0;
	    alpha.elem(1,1) = 1.0;
	    alpha.elem(2,2) = 1.0 + f;
	    alpha.elem(0,2) = f;
	  }
	  else {

	  }
	}
	else if (isym==12){
	  if (csystem_sub==0){
	    /* Result: E - E0 = 0.5 * V0 * (4.0 * C44 + 4.0 * C66 + 8.0 * C46) * f^2 */
	    alpha.elem(0,0) = 1.0;
	    alpha.elem(1,1) = 1.0;
	    alpha.elem(2,2) = 1.0;
	    alpha.elem(1,2) = f;
	    alpha.elem(0,1) = f;
	  }
	  else {


	  }
	}


	else quit=true;
      }
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      else if (csystem=="trigonal"){



      }
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      else if (csystem=="tetragonal"){



      }
      // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      else quit=true;
#endif



#if 0


  
  if (csystem=="cubic"){
    C44 = Clincomb.elem(4,4) / 2.0;
    C11 = Clincomb.elem(1,1) * 2.0;
    C12 = Clincomb.elem(1,2);
  }
  else if (csystem=="hexagonal"){


    C33 = 2.0 * Clincomb[2];
    C11 = (Clincomb[0] + Clincomb[1])/2.0;
    C12 = Clincomb[0] - C11;
    C13 = (2.0*Clincomb[3] - C11 - C33)/2.0;
    C44 = Clincomb[4]/2.0;
    C.elem(0,0) = C11;
    C.elem(0,1) = C12;
    C.elem(3,3) = C44;
    C.elem(2,2) = C33;
    C.elem(0,2) = C13;
  }
  else if (csystem=="orthorombic"){
    C11 = 2.0 * Clincomb[0];
    C22 = 2.0 * Clincomb[1];
    C33 = 2.0 * Clincomb[2];
    C12 = 0.5 * ( 2.0 * Clincomb[3] - C11 - C22 );
    C13 = 0.5 * ( 2.0 * Clincomb[4] - C11 - C33 );
    C23 = 0.5 * ( 2.0 * Clincomb[5] - C22 - C33 );
    C44 = Clincomb[6]/2.0;
    C55 = Clincomb[7]/2.0;
    C66 = Clincomb[8]/2.0;
    C.elem(0,0) = C11;
    C.elem(1,1) = C22;
    C.elem(2,2) = C33;
    C.elem(0,1) = C12;
    C.elem(0,2) = C13;
    C.elem(1,2) = C23;
    C.elem(3,3) = C44;
    C.elem(4,4) = C55;
    C.elem(5,5) = C66;
  }
  else if (csystem=="triclinic"){
    C11 = 2.0 * Clincomb[0];
    C22 = 2.0 * Clincomb[1];
    C33 = 2.0 * Clincomb[2];
    C44 = 4.0 * Clincomb[3];
    C55 = 4.0 * Clincomb[4];
    C66 = 4.0 * Clincomb[5];

    C12 = (2.0 * Clincomb[6]  - C11 - C22)/2.0;
    C13 = (2.0 * Clincomb[7]  - C11 - C33)/2.0;
    C14 = (2.0 * Clincomb[8]  - C11 - 4.0 * C44)/2.0;
    C15 = (2.0 * Clincomb[9]  - C11 - 4.0 * C55)/2.0;
    C16 = (2.0 * Clincomb[10] - C11 - 4.0 * C66)/2.0;

    C23 = (2.0 * Clincomb[11]  - C22 - C33)/2.0;
    C24 = (2.0 * Clincomb[12]  - C22 - 4.0 * C44)/2.0;
    C25 = (2.0 * Clincomb[13]  - C22 - 4.0 * C55)/2.0;
    C26 = (2.0 * Clincomb[14]  - C22 - 4.0 * C66)/2.0;

    C34 = (2.0 * Clincomb[15]  - C33 - 4.0 * C44)/2.0;
    C35 = (2.0 * Clincomb[16]  - C33 - 4.0 * C55)/2.0;
    C36 = (2.0 * Clincomb[17]  - C33 - 4.0 * C66)/2.0;

    C45 = (2.0 * Clincomb[18]  - 4.0 * C44 - 4.0 * C55)/8.0;
    C46 = (2.0 * Clincomb[19]  - 4.0 * C44 - 4.0 * C66)/8.0;

    C56 = (2.0 * Clincomb[20]  - 4.0 * C55 - 4.0 * C66)/8.0;

    C.elem(0,0) = C11;
    C.elem(1,1) = C22;
    C.elem(2,2) = C33;
    C.elem(3,3) = C44;
    C.elem(4,4) = C55;
    C.elem(5,5) = C66;

    C.elem(0,1) = C12;
    C.elem(0,2) = C13;
    C.elem(0,3) = C14;
    C.elem(0,4) = C15;
    C.elem(0,5) = C16;

    C.elem(1,2) = C23;
    C.elem(1,3) = C24;
    C.elem(1,4) = C25;
    C.elem(1,5) = C26;

    C.elem(2,3) = C34;
    C.elem(2,4) = C35;
    C.elem(2,5) = C36;

    C.elem(3,4) = C45;
    C.elem(3,5) = C46;

    C.elem(4,5) = C56;
  }
  else if (csystem=="monoclinic"){
    C11 = 2.0 * Clincomb[0];
    C22 = 2.0 * Clincomb[1];
    C33 = 2.0 * Clincomb[2];
    C44 = 4.0 * Clincomb[3];
    C55 = 4.0 * Clincomb[4];
    C66 = 4.0 * Clincomb[5];

    C12 = (2.0 * Clincomb[6]  - C11 - C22)/2.0;
    C13 = (2.0 * Clincomb[7]  - C11 - C33)/2.0;

    C15 = (2.0 * Clincomb[8]  - C11 - 4.0 * C55)/2.0;

    C23 = (2.0 * Clincomb[9]  - C22 - C33)/2.0;
    C25 = (2.0 * Clincomb[10]  - C22 - 4.0 * C55)/2.0;

    C35 = (2.0 * Clincomb[11]  - C33 - 4.0 * C55)/2.0;

    C46 = (2.0 * Clincomb[12]  - 4.0 * C44 - 4.0 * C66)/8.0;

    C.elem(0,0) = C11;
    C.elem(1,1) = C22;
    C.elem(2,2) = C33;
    C.elem(3,3) = C44;
    C.elem(4,4) = C55;
    C.elem(5,5) = C66;

    C.elem(0,1) = C12;
    C.elem(0,2) = C13;

    C.elem(0,4) = C15;

    C.elem(1,2) = C23;

    C.elem(1,4) = C25;

    C.elem(2,4) = C35;

    C.elem(3,5) = C46;
  }



  C.elem(0,0) = C11;
  C.elem(0,1) = C12;
  C.elem(0,2) = C13;
  C.elem(0,3) = C14;
  C.elem(0,4) = C15;
  C.elem(0,5) = C16;

  C.elem(1,1) = C22;
  C.elem(1,2) = C23;
  C.elem(1,3) = C24;
  C.elem(1,4) = C25;
  C.elem(1,5) = C26;

  C.elem(2,2) = C33;
  C.elem(2,3) = C34;
  C.elem(2,4) = C35;
  C.elem(2,5) = C36;

  C.elem(3,3) = C44;
  C.elem(3,4) = C45;
  C.elem(3,5) = C46;

  C.elem(4,4) = C55;
  C.elem(4,5) = C56;

  C.elem(5,5) = C66;


#endif



