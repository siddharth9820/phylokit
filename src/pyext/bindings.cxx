#include <boost/python.hpp>
#include <boost/python/numeric.hpp>
#include "../newick.hpp"
#include "../TaxonSet.hpp"
#include <numpy/arrayobject.h>


using namespace boost::python;

TaxonSet* newick_to_taxon_set(boost::python::list& strings) {
  
  unordered_set<string> taxa;

  for (int i = 0; i < len(strings); i++) {
    newick_to_ts(extract<string>(strings[i]), taxa);
  }

  TaxonSet* ts = new TaxonSet(taxa.size());

  for (string s : taxa) {
    ts->add(s);
  }
  ts->freeze();
  
  return ts;
}

void newick_to_dm_(const string& tree, TaxonSet& ts , object& dists, object& masks) {

  dm_type dist_mat(boost::extents[ts.size()][ts.size()]);
  dm_type mask_mat(boost::extents[ts.size()][ts.size()]);
  
  newick_to_dm(tree, ts, dist_mat, mask_mat);
  
  for (int i = 0; i < ts.size(); i++) {
    for (int j = 0; j < ts.size(); j++) {
      *(double*)PyArray_GETPTR2(dists.ptr(), i, j) = dist_mat[i][j];
      *(double*)PyArray_GETPTR2(masks.ptr(), i, j) = mask_mat[i][j]; 
    }
  }
  
}

boost::python::list newick_to_tax_indices_(const string& tree, TaxonSet& ts) {
  boost::python::list indices;
  unordered_set<string> taxa;

  newick_to_ts(tree, taxa);
  
  for (const string& t : taxa) {
    indices.append(str(ts[t]));
  }    
  
  return indices;
}

BOOST_PYTHON_MODULE(libphylokit)
{

  class_<TaxonSet>("TaxonSet", init<int>())
    .def("size", &TaxonSet::size)
    .def("__getitem__", &TaxonSet::get, return_value_policy<copy_const_reference>());
  
  class_<dm_type>("distmat");
  
  def("newick_to_taxon_set", newick_to_taxon_set, return_value_policy<manage_new_object>());
  def("newick_to_tax_indices", newick_to_tax_indices_);

  def("newick_to_dm", newick_to_dm_);
  
}
