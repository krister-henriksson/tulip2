

void MDSystem::handle_pbc_of_positions(const double lowlim){
  Vector<double> drs(3), drc(3);
  int i;

  double pf1=-0.5, pf2=0.5;
  if (lowlim>0){
    pf1 = 0.0;
    pf2 = 1.0;
  }

  for (i=0; i<natoms(); ++i){

    drc[0] = pos[i][0];
    drc[1] = pos[i][1];
    drc[2] = pos[i][2];

    // Get distance in skew coordinate system, where periodics can be checked:
    drs[0] = drc[0];
    drs[1] = drc[1];
    drs[2] = drc[2];
    if (! isCart){
      drs[0] = Bravaismatrix_inv.elem(0,0) * drc[0]
	+ Bravaismatrix_inv.elem(0,1) * drc[1]
	+ Bravaismatrix_inv.elem(0,2) * drc[2];
      drs[1] = Bravaismatrix_inv.elem(1,0) * drc[0]
	+ Bravaismatrix_inv.elem(1,1) * drc[1]
	+ Bravaismatrix_inv.elem(1,2) * drc[2];
      drs[2] = Bravaismatrix_inv.elem(2,0) * drc[0]
	+ Bravaismatrix_inv.elem(2,1) * drc[1]
	+ Bravaismatrix_inv.elem(2,2) * drc[2];
    }
      
    // Periodics check:
    while (pbc[0] && drs[0] <  pf1 * boxlen[0]) drs[0] += boxlen[0];
    while (pbc[0] && drs[0] >= pf2 * boxlen[0]) drs[0] -= boxlen[0];
    while (pbc[1] && drs[1] <  pf1 * boxlen[1]) drs[1] += boxlen[1];
    while (pbc[1] && drs[1] >= pf2 * boxlen[1]) drs[1] -= boxlen[1];
    while (pbc[2] && drs[2] <  pf1 * boxlen[2]) drs[2] += boxlen[2];
    while (pbc[2] && drs[2] >= pf2 * boxlen[2]) drs[2] -= boxlen[2];

    // Get distance in Cartesian coordinate system:
    drc[0] = drs[0];
    drc[1] = drs[1];
    drc[2] = drs[2];
    if (! isCart){
      drc[0] = boxdir.elem(0,0) * drs[0] + boxdir.elem(0,1) * drs[1] + boxdir.elem(0,2) * drs[2];
      drc[1] = boxdir.elem(1,0) * drs[0] + boxdir.elem(1,1) * drs[1] + boxdir.elem(1,2) * drs[2];
      drc[2] = boxdir.elem(2,0) * drs[0] + boxdir.elem(2,1) * drs[1] + boxdir.elem(2,2) * drs[2];
    }
    pos[i][0] = drc[0];
    pos[i][1] = drc[1];
    pos[i][2] = drc[2];
      
  }
}



