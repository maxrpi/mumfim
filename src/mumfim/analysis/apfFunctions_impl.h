namespace amsi
{
  template <typename I>
    void setEntitiesNode(apf::Field * fld, double vl, I bgn, I nd)
  {
    for(auto ent = bgn; ent != nd; ++ent)
    {
      apf::MeshEntity * mnt = reinterpret_cast<apf::MeshEntity*>(*ent);
      apf::setScalar(fld,mnt,0,vl);
    }
  }
}
