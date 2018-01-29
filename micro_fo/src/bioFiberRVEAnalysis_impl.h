namespace bio
{
  template <typename I>
    void applyRVEBC(I bnd_bgn, I bnd_end, apf::Numbering * nm)
  {
    for(auto ent = bnd_bgn; ent != bnd_end; ++ent)
    {
      for(int ii = 0; ii < 3; ++ii)
        apf::fix(nm,*ent,0,ii,true);
      /*
      apf::NewArray<int> dofs;
      int nedofs = apf::getElementNumbers(nm,*ent,dofs);
      apf::DynamicVector zeroes(nedofs);
      zeroes.zero();
      ops->set(f,nedofs,&dofs[0],&zeroes[0]);
      apf::DynamicMatrix eye(nedofs,nedofs);
      eye.zero();
      for(int ii = 0; ii < nedofs; ++ii)
        for(int jj = 0; jj < nedofs; ++jj)
          eye(ii,jj) = 1.0;
      ops->set(k,nedofs,&dofs[0],nedofs,&dofs[0],&eye(0,0));
      */
    }
  }
}
