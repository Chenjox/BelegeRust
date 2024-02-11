use nalgebra::{constraint, Dyn, OMatrix, Point2, Vector2, Vector3};

use crate::svd_helper::get_nullspace;


mod svd_helper;

fn main(){
  println!("Hello World");

  let lines = vec![Point2::new(0.0,1.0)];
  let points = vec![Point2::new(0.,0.), Point2::new(1.,0.0),Point2::new(2.,0.)];
  let constrained = vec![false,false,false,true];

  // Alles in homogenen Coordinaten
  let mut hom_points: Vec<Vector3<f64>> = points.iter().map(|p| p.to_homogeneous().normalize()).collect();
  let mut hom_lines: Vec<Vector3<f64>> = lines.iter().map(|p| {
    let mut new_p = p.to_homogeneous();
    new_p[2] = 0.0;
    new_p.normalize()
    }
  ).collect();

  hom_points.append(&mut hom_lines);

  let hom_points = hom_points;

  let edges = vec![[0,1],[1,2],[2,3],[0,3]];
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
  if let Some((left, right)) = get_nullspace(&reduced_rigidity_matrix) {

    let null = right.transpose();
    let null_space = null.rows(rank, null.nrows()-rank);

    //println!("{:3.3},{:3.3}",left,null);
    //println!("{:3.3}",null_space);
    let null_space_size = null.nrows() - rank;

    let mut solution_matrix_real_space = OMatrix::<f64,Dyn,Dyn>::zeros(null_space_size, dim*(num_vertices-num_constrained));

    for point_index in 0..(num_vertices-num_constrained) {
      let view = null_space.columns_with_step(point_index*projective_dim, projective_dim, 0);

      let point = hom_points[point_index];
      let point = point / point[2];

      let velo = view*point.norm();
      println!("{:3.3}",velo);
      for i in 0..null_space_size {
        for j in 0..dim {
          solution_matrix_real_space[(i,point_index*dim+j)] = velo[(i,j)];
        }
      }
    }

    // Das sind alle Basisvektoren der Mechanismen des Tragwerks
    println!("{:3.3}",solution_matrix_real_space);
    // die sind linear unabhängig, aber nicht zwangsweise orthogonal
    // ergo: basiswechsel ist angesagt.
    let mut angle_matrix = OMatrix::<f64,Dyn,Dyn>::zeros(null_space_size,null_space_size);
    for i in 0..null_space_size {
      for j in 0..null_space_size {
        angle_matrix[(i,j)] = solution_matrix_real_space.row(i).dot(&solution_matrix_real_space.row(j));
      }
    }
    // Diese Matrix wird _immer_ positiv semidefinit sein und symmetrisch = LU oder Cholesky
    println!("{:3.3}",angle_matrix);
    let chol = angle_matrix.symmetric_eigen();
    let orthogonalest_mechanism = chol.eigenvectors*solution_matrix_real_space;
    println!("{:3.3}",orthogonalest_mechanism)
  };

}
