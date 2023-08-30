
fn get_tragwerk() -> Fachwerk2D {
  let mut v = Vec::new();

  // Knoten (noch nicht in Koordinaten die ich brauch)
  for i in 0..7 {
    v.push(Point2D::new_unconstrained(i + 1, 3.0 * i as f64, 7.0));
  }
  for i in 0..7 {
    if (10..=12).contains(&(i + 8)) {
      v.push(Point2D::new_unconstrained(i + 8, 3.0 * i as f64, 4.5));
    } else {
      v.push(Point2D::new_unconstrained(i + 8, 3.0 * i as f64, 3.0));
    }
  }
  for i in 0..2 {
    v.push(Point2D::new(
      i + 15,
      3.0 * i as f64,
      0.0,
      vec![CONSTRAINS::XRestrain, CONSTRAINS::YRestrain],
    ));
    v.push(Point2D::new(
      18 - i,
      3.0 * (1 - i) as f64 + 15.0,
      0.0,
      vec![CONSTRAINS::XRestrain, CONSTRAINS::YRestrain],
    ));
  }
  // Kanten
  let mut kant = Vec::new();
  for i in 0..7 {
    kant.push(Beam2D::new(i + 1, i + 8, 1.0, 1.0));
  }
  for i in 0..6 {
    kant.push(Beam2D::new(i + 1, i + 2, 1.0, 1.0));
    kant.push(Beam2D::new(i + 8, i + 9, 1.0, 1.0));
  }
  for i in 0..3 {
    kant.push(Beam2D::new(i + 1, i + 9, 1.0, 1.0));
    kant.push(Beam2D::new(7 - i, 13 - i, 1.0, 1.0));
  }
  for i in 0..2 {
    kant.push(Beam2D::new(i + 8, i + 15, 1.0, 1.0));
    kant.push(Beam2D::new(14 - i, 18 - i, 1.0, 1.0));
  }
  kant.push(Beam2D::new(8, 16, 1.0, 1.0));
  kant.push(Beam2D::new(14, 17, 1.0, 1.0));

  return Fachwerk2D::new(v, kant);
}

fn get_loading() -> Vec<Load2D> {
  return vec![
    Load2D::new_from_components(1, 3.0, 0.),
    Load2D::new_from_components(4, 0., -2.0),
    Load2D::new_from_components(5, 0., -1.0),
  ];
}

fn get_testtragwerk() -> Fachwerk2D {
  let points = vec![
    Point2D::new(
      0,
      0.,
      0.,
      vec![CONSTRAINS::XRestrain, CONSTRAINS::YRestrain],
    ),
    Point2D::new_unconstrained(1, 1., 0.),
    Point2D::new_unconstrained(2, 2., 0.),
    Point2D::new(
      3,
      3.,
      0.,
      vec![CONSTRAINS::XRestrain, CONSTRAINS::YRestrain],
    ),
  ];
  let beams = vec![
    Beam2D::new(0, 1, 10., 1.0),
    Beam2D::new(1, 2, 10., 1.0),
    Beam2D::new(2, 3, 10., 1.0),
  ];

  return Fachwerk2D::new(points, beams);
}

fn get_loading_test2() -> Vec<Load2D> {
  return vec![
    Load2D::new_from_components(1, 3.0, 0.),
    Load2D::new_from_components(4, 0., -2.0),
    Load2D::new_from_components(5, 0., -1.0),
  ];
}

fn get_testtragwerk2() -> Fachwerk2D {
  let points = vec![
    Point2D::new(
      1,
      0.,
      0.,
      vec![CONSTRAINS::XRestrain, CONSTRAINS::YRestrain],
    ),
    Point2D::new(2, 2., 0., vec![CONSTRAINS::YRestrain]),
    Point2D::new_unconstrained(3, 1., 2.),
    Point2D::new_unconstrained(4, 1., 0.),
    Point2D::new_unconstrained(5, 0.5, 1.),
    Point2D::new_unconstrained(6, 1.5, 1.),
  ];
  let beams = vec![
    Beam2D::new(1, 4, 10., 1.0),
    Beam2D::new(1, 5, 10., 1.0),
    Beam2D::new(2, 4, 10., 1.0),
    Beam2D::new(2, 6, 10., 1.0),
    Beam2D::new(3, 4, 10., 1.0),
    Beam2D::new(3, 5, 10., 1.0),
    Beam2D::new(3, 6, 10., 1.0),
    Beam2D::new(4, 5, 10., 1.0),
    Beam2D::new(4, 6, 10., 1.0)
  ];

  return Fachwerk2D::new(points, beams);
}