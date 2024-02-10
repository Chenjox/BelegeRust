use nalgebra::{constraint, Dyn, OMatrix, Point2, Vector3};




fn main(){
  println!("Hello World");


  let points = vec![Point2::new(0.,0.), Point2::new(1.,0.0), Point2::new(2.,0.0),Point2::new(3.,0.)];
  let constrained = vec![true,false,false,true];

  // Alles in homogenen Coordinaten
  let hom_points: Vec<Vector3<f64>> = points.iter().map(|p| p.to_homogeneous().normalize()).collect();

  let edges = vec![[0,1],[1,2],[2,3]];
  let dim = 2;

  let num_vertices = hom_points.len();
  let num_edges = edges.len();
  let projective_dim = dim+1;

  let num_constrained: usize = constrained.iter().map(|f| if *f { 1 } else { 0 }).sum();

  let mut rigidity_matrix = OMatrix::<f64,Dyn,Dyn>::zeros(num_edges+num_vertices, projective_dim*num_vertices);

  for p in &hom_points {
    println!("{}",p)
  }

  // ZEILEN Vektoren, nicht spalten!
  for (index,edge) in edges.iter().enumerate() {
    let from_point = &hom_points[edge[0]];
    let to_point = &hom_points[edge[1]];
    
    for i in 0..projective_dim {
      rigidity_matrix[(index,edge[1]*projective_dim+i)] = to_point[i] - from_point[i];
      rigidity_matrix[(index,edge[0]*projective_dim+i)] = -(to_point[i] - from_point[i]);
    }
  }
  
  // Punkte einfügen
  for (index,point) in hom_points.iter().enumerate() {
    for i in 0..projective_dim {
      rigidity_matrix[(num_edges+index,index*projective_dim+i)] = point[i];
    }
  }

  // jetzt bauen wir die constraints ein:
  // und nehmen an, dass alle Edges dazugehören
  let mut constraint_points_seen = 0;
  let mut reduced_rigidity_matrix =  OMatrix::<f64,Dyn,Dyn>::zeros(num_edges+(num_vertices-num_constrained), projective_dim*(num_vertices-num_constrained));
  for i in 0..hom_points.len() {
    if !constrained[i] {
      for d in 0..projective_dim {
        // Kantenconstraints
        for edge_index in 0..num_edges {
          reduced_rigidity_matrix[(edge_index,constraint_points_seen*projective_dim+d)] = rigidity_matrix[(edge_index,i*projective_dim+d)]
        }
        // Punktwert selbst
        reduced_rigidity_matrix[(num_edges+constraint_points_seen,constraint_points_seen*projective_dim+d)] = rigidity_matrix[(num_edges+i,i*projective_dim+d)]
        
      }
      constraint_points_seen += 1;
    }
  }

  println!("{:3.3}",rigidity_matrix);
  println!("{}",reduced_rigidity_matrix);

  // irgendwie kommt der nullraum nicht raus
  reduced_rigidity_matrix.iter_mut().for_each(|f| if f.abs() < 1e-10 { *f = 0.0 });

  println!("R = {:3.3}",reduced_rigidity_matrix);


  let rank = reduced_rigidity_matrix.rank(1e-10);
  let svd = reduced_rigidity_matrix.transpose().svd(true, true);
  let qr = reduced_rigidity_matrix.transpose().qr();
  let null = qr.q().transpose();
  let rot = qr.r();
  //let sing = svd.singular_values;
  //let null = svd.u.unwrap();


  //println!("{:3.3}",sing);
  println!("Q = {:3.3}",null);
  println!("R = {:3.3}",rot);

  let null_space = null.rows(rank, null.nrows()-rank);

  println!("{:3.3}",null_space);
}