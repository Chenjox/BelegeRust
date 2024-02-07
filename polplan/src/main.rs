use nalgebra::{Dyn, OMatrix, Point2, Vector3};




fn main(){
  println!("Hello World");


  let points = vec![Point2::new(0.,0.), Point2::new(1.,0.), Point2::new(2.,0.)];

  // Alles in homogenen Coordinaten
  let hom_points: Vec<Vector3<f64>> = points.iter().map(|p| p.to_homogeneous().normalize()).collect();

  let edges = vec![[0,1],[1,2],[0,2]];
  let dim = 2;

  let num_vertices = hom_points.len();
  let num_edges = edges.len();
  let projective_dim = dim+1;


  let mut rigidity_matrix = OMatrix::<f64,Dyn,Dyn>::zeros(num_edges+num_vertices, projective_dim*num_vertices);

  for p in &hom_points {
    println!("{}",p)
  }

  for (index,edge) in edges.iter().enumerate() {
    let from_point = &hom_points[edge[0]];
    let to_point = &hom_points[edge[1]];

    for i in 0..projective_dim {
      rigidity_matrix[(index,edge[0]*projective_dim+i)] = to_point[i];
      rigidity_matrix[(index,edge[1]*projective_dim+i)] = from_point[i];
    }
  }

  // Punkte einfügen
  for (index,point) in hom_points.iter().enumerate() {
    for i in 0..projective_dim {
      rigidity_matrix[(num_edges+index,index*projective_dim+i)] = point[i];
    }
  }

  println!("{:3.3}",rigidity_matrix);

  let rank = rigidity_matrix.rank(1e-10);
  let svd = rigidity_matrix.svd(true, true);

  let mechanism = svd.v_t.unwrap();

  let sing = svd.singular_values;
  let num_singular = sing.len();
  println!("{},{}",num_singular,rank);
  println!("{:3.3}",mechanism);

  let null_space = mechanism.rows(rank, num_singular-rank);

  println!("{:3.3}",null_space);

  for point_index in 0..num_vertices {
    let view = null_space.columns_with_step(point_index*projective_dim, projective_dim, 0);

    let view = &view.transpose();
    let projection = Vector3::new(0.0,0.0,1.);
    let point = hom_points[point_index];

    let p = (view - view.dot(&point) * projection)/point.norm();
    println!("{:3.3}",p);
  }

}
