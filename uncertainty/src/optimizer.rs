use rand::prelude::*;

fn add(a: f64, av: &Vec<f64>, b: f64, bv: &Vec<f64>, erg: &mut Vec<f64>) {
  for i in 0..erg.len() {
    erg[i] = a * av[i] + b * bv[i]
  }
}

fn simplex_centroid(vals: &[Vec<f64>], excluded: usize, n: usize) -> Vec<f64> {
  let mut g = vec![0.0; n];
  for j in 0..n {
    for i in 0..n + 1 {
      if i != excluded {
        g[j] = g[j] + vals[i][j];
      }
    }
    g[j] = 1.0 / ((n) as f64) * g[j];
  }
  return g;
}

fn is_in_bounds(
  lower_bound: &Vec<f64>,
  upper_bound: &Vec<f64>,
  value: &Vec<f64>,
  n: usize,
) -> bool {
  let mut in_bound = true;
  for j in 0..n {
    in_bound = in_bound && lower_bound[j] <= value[j] && upper_bound[j] >= value[j];
  }
  return in_bound;
}

fn reflection(center: &Vec<f64>, tobereflected: &Vec<f64>, n: usize) -> Vec<f64> {
  let mut g = vec![0.0; n];
  for j in 0..n {
    g[j] = center[j] + 1.0 * (center[j] - tobereflected[j]);
  }
  return g;
}

fn expand(center: &Vec<f64>, tobe: &Vec<f64>, n: usize) -> Vec<f64> {
  let mut g = vec![0.0; n];
  for j in 0..n {
    g[j] = center[j] + 2.0 * (tobe[j] - center[j]);
  }
  return g;
}

fn contract(center: &Vec<f64>, tobe: &Vec<f64>, n: usize) -> Vec<f64> {
  let mut g = vec![0.0; n];
  for j in 0..n {
    g[j] = center[j] + 0.5 * (tobe[j] - center[j]);
  }
  return g;
}

fn sort_vals(points: &mut Vec<Vec<f64>>, point_values: &mut Vec<f64>) {
  let mut i = 1;
  while i < point_values.len() {
    let x = point_values[i];
    let px = points[i].clone();
    let mut j = i - 1;
    let mut zero_signal = false;
    'inner: while point_values[j] > x {
      if j == 0 {
        point_values[j + 1] = point_values[j];
        points[j + 1] = points[j].clone();
        zero_signal = true;
        break 'inner;
      } else {
        point_values[j + 1] = point_values[j];
        points[j + 1] = points[j].clone();
        j = j - 1;
      }
    }
    if zero_signal {
      point_values[0] = x;
      points[0] = px;
    } else {
      point_values[j + 1] = x;
      points[j + 1] = px;
    }

    i = i + 1;
  }
}

pub fn simplex_optimization<F>(
  lower_bound: &Vec<f64>,
  upper_bound: &Vec<f64>,
  start_point: &Vec<f64>,
  func: F,
) -> Vec<f64>
where
  F: Fn(&Vec<f64>) -> f64,
{
  let n = lower_bound.len();

  let mut points = Vec::with_capacity(n + 1);
  // Anfangspunkte
  {
    // Werden mittels lerp erstellt.
    for i in 0..n {
      // nehme einen Schritt im Bereich
      let interval = (start_point[i] - lower_bound[i]) * 0.5;
      // erstelle einen neuen punkt
      let mut new_point = vec![0.0; n];
      for k in 0..n {
        // der in der Koordinate k vom startpunkt abweicht
        if k == i {
          new_point[k] = start_point[k] + interval * 0.7
        } else {
          // und sonst gleich ist.
          new_point[k] = start_point[k]
        }
      }
      // und füge Ihn zur der point liste hinzu
      points.push(new_point);
    }
  }
  // letztlich der startpunkt
  points.push(start_point.clone());
  // erste funktionswerte
  let mut point_values = vec![0.0; n + 1];
  for i in 0..n + 1 {
    point_values[i] = func(&points[i]);
  }
  // sortieren der Punkte
  let mut iter_counter = 0;
  let MAX_ITER = 10000;
  'itera: loop {
    if iter_counter > MAX_ITER {
      break 'itera;
    }
    iter_counter += 1;
    sort_vals(&mut points, &mut point_values);
    // hier kommt die terminierung hin.
    {
      let mut dist = 0.0;
      for k in 0..n {
        dist += (points[0][k] - points[n][k]).powi(2);
      }
      dist = dist.sqrt();
      if dist < 1e-10 {
        //println!("{}", iter_counter);
        break 'itera;
      }
    }

    // vorne steht das beste, hinten das schlimmste
    // Schritt 2
    let o = simplex_centroid(&points, n, n);

    // Schritt 3
    let r = reflection(&o, &points[n], n);
    let fr = if is_in_bounds(&lower_bound, &upper_bound, &r, n) {
      func(&r)
    } else {
      point_values[n] + 1.0
    };

    if point_values[0] <= fr && fr < point_values[n - 1] {
      points[n] = r;
      point_values[n] = fr;
      continue 'itera;
    }
    // Schritt 4.
    if fr < point_values[0] {
      // reflected is best!
      let e = expand(&o, &r, n);
      let fe = if is_in_bounds(&lower_bound, &upper_bound, &e, n) {
        func(&e)
      } else {
        point_values[n] + 1.0
      };
      if fe < fr {
        points[n] = e;
        point_values[n] = fe;
      } else {
        points[n] = r;
        point_values[n] = fr;
      }
      continue 'itera;
    }
    // Schritt 5
    if fr < point_values[n] {
      let c = contract(&o, &r, n);
      let fc = if is_in_bounds(&lower_bound, &upper_bound, &c, n) {
        func(&c)
      } else {
        point_values[n] + 1.0
      };
      if fc < fr {
        points[n] = c;
        point_values[n] = fc;
        continue 'itera;
      } else {
        // schrinking
        for i in 1..n + 1 {
          for k in 0..n {
            points[i][k] = points[0][k] + 0.5 * (points[i][k] - points[0][k]);
          }
        }
        for i in 1..n + 1 {
          point_values[i] = func(&points[i]);
        }
        continue 'itera;
      }
    }
    if fr >= point_values[n] {
      let c = contract(&o, &points[n], n); // Dieser Punkt muss innerhalb liegen!
      let fc = func(&c);
      if fc < point_values[n] {
        points[n] = c;
        point_values[n] = fc;
        continue 'itera;
      } else {
        // schrinking
        for i in 1..n + 1 {
          for k in 0..n {
            points[i][k] = points[0][k] + 0.5 * (points[i][k] - points[0][k]);
          }
        }
        for i in 1..n + 1 {
          point_values[i] = func(&points[i]);
        }
        continue 'itera;
      }
    }
  }
  return points[0].clone();
}

pub fn such_optimierung<F>(
  lower_bound: &Vec<f64>,
  upper_bound: &Vec<f64>,
  start_point: &Vec<f64>,
  func: F,
) -> Vec<f64>
where
  F: Fn(&Vec<f64>) -> f64,
{
  let mut rng = rand::thread_rng();
  let anzahl_kinder = 20;

  let dimension = start_point.len();

  // Schrittgrößen für jede Dimension
  let mut step_size = vec![0.; dimension];
  for i in 0..dimension {
    step_size[i] = (upper_bound[i] - lower_bound[i]) * 0.5;
  }
  let step_size = step_size;

  // Schrittfaktor für jede Dimension
  let mut step_factor = 1.0;
  //

  let mut current_point = start_point.clone();
  let mut current_best = func(&start_point);

  loop {
    // Erstellen der unterschiedlichen Richtungen
    let mut iter_result = Vec::new();
    for i in 0..anzahl_kinder {
      let mut direction = vec![0.; dimension]; // Richtungsvektor
      let mut norm = 0.; // Normierungssumme
      for d in 0..dimension {
        let coord: f64 = rng.gen();
        direction[d] = coord;
        norm += coord * coord;
      }
      norm = norm.sqrt(); // Hey look, es ist eine Norm!
      for d in 0..dimension {
        // Normierung, Größenskalierung und verkleinerungsfaktor in einem!
        direction[d] = direction[d] / norm * step_size[d] * step_factor;
      }
      // Jetzt der neue Punkt
      let mut new_point = vec![0.; dimension];
      for d in 0..dimension {
        new_point[d] = current_point[d] + direction[d];
      }
      //Randprojektion
      project_onto_domain(lower_bound, upper_bound, &mut new_point);

      let function_value = func(&new_point);
      // Abspeichern der Kinder
      iter_result.push((new_point,function_value));
    }
    let mut better_count = 0;
    // Durchzählen
    for i in &iter_result {
      if i.1 <= current_best {
        better_count += 1;
      }
    }
    // Erfolg, oder nicht erfolg
    if better_count as f64/anzahl_kinder as f64 >= 0.1 {
      iter_result.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
      let (point, value) = &iter_result[0];
      current_point = point.clone();
      // Abbruchbedingung:
      if (current_best - value).abs() < 1e-15 {
        break;
      }
      current_best = *value;
    } else {
      step_factor *= 0.5;
    }
  }

  return current_point;
}

fn project_onto_domain(lower_bound: &Vec<f64>, upper_bound: &Vec<f64>, point: &mut Vec<f64>) {
  let n = point.len();
  for i in 0..n {
    if lower_bound[i] > point[i] {
      point[i] = lower_bound[i]
    }
    if upper_bound[i] < point[i] {
      point[i] = upper_bound[i]
    }
  }
}
