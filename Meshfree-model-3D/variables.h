dimension pEBC(3,1000),npEBC(4,1000),npNBC(4,1000),pNBC(3,1000)
dimension Dmat(6,6)
dimension x(nx,numd),noCell1(ng,ncn),ds(nx,numd)
dimension xc(nx,numdq),noCell(ng,numc)
dimension gauss(2,nqc),gs(ng,numg)
dimension gpos(nx),nv(numd),ph(10,numd)
dimension ak(3*numd,3*numd),GSPk(9*numd*numd)
dimension ne(3*numd),force(3*numd)
dimension u2(nx,numd),disp(3*numd)
dimension Stressnode(6,numd)
common/para/xlength,ylength,zlength,p,young,anu,aimo
common/rpim/ALFC,DC,Q,nRBF
common /basis/mbasis
