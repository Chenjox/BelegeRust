use std::f64::consts;

#[derive(Clone)]
pub struct Schraubenfestigkeitsklasse {
  klasse: [u32; 2],
}
impl Schraubenfestigkeitsklasse {
  pub fn get_festigkeitsklassen_deutschland() -> [[u32; 2]; 4] {
    [[4, 6], [5, 6], [8, 8], [10, 9]]
  }
  /// In N/mm^2
  pub fn zugfestigkeit(&self) -> f64 {
    return self.klasse[0] as f64 * 100.0;
  }
  /// In N/mm^2
  pub fn streck_dehn_grenze(&self) -> f64 {
    return (self.klasse[0] * self.klasse[1] * 10) as f64;
  }
  pub fn is_streck_oder_dehngrenz(&self) -> bool {
    return self.klasse[0] >= 8;
  }
}

pub type SFK = Schraubenfestigkeitsklasse;

#[derive(Clone)]
pub struct ISOSchraube {
  festigkeitsklasse: SFK,
  schraubengroesse: f64,
  gewindesteigung: f64,
  passschraube: bool,
}

impl ISOSchraube {
  pub fn new(sfk: [u32; 2], schraubengroesse: u32, passschraube: bool) -> Self {
    match schraubengroesse {
      1 => ISOSchraube {
        festigkeitsklasse: SFK { klasse: sfk },
        schraubengroesse: 1.0,
        gewindesteigung: 0.25,
        passschraube,
      },
      2 => ISOSchraube {
        festigkeitsklasse: SFK { klasse: sfk },
        schraubengroesse: 2.0,
        gewindesteigung: 0.4,
        passschraube,
      },
      3 => ISOSchraube {
        festigkeitsklasse: SFK { klasse: sfk },
        schraubengroesse: 3.0,
        gewindesteigung: 0.5,
        passschraube,
      },
      12 => ISOSchraube {
        festigkeitsklasse: SFK { klasse: sfk },
        schraubengroesse: 12.0,
        gewindesteigung: 1.75,
        passschraube,
      },
      16 => ISOSchraube {
        festigkeitsklasse: SFK { klasse: sfk },
        schraubengroesse: 16.0,
        gewindesteigung: 2.,
        passschraube,
      },
      20 => ISOSchraube {
        festigkeitsklasse: SFK { klasse: sfk },
        schraubengroesse: 20.0,
        gewindesteigung: 2.5,
        passschraube,
      },
      22 => ISOSchraube {
        festigkeitsklasse: SFK { klasse: sfk },
        schraubengroesse: 22.0,
        gewindesteigung: 2.5,
        passschraube,
      },
      24 => ISOSchraube {
        festigkeitsklasse: SFK { klasse: sfk },
        schraubengroesse: 24.0,
        gewindesteigung: 3.,
        passschraube,
      },
      27 => ISOSchraube {
        festigkeitsklasse: SFK { klasse: sfk },
        schraubengroesse: 27.0,
        gewindesteigung: 3.,
        passschraube,
      },
      30 => ISOSchraube {
        festigkeitsklasse: SFK { klasse: sfk },
        schraubengroesse: 30.0,
        gewindesteigung: 3.5,
        passschraube,
      },
      36 => ISOSchraube {
        festigkeitsklasse: SFK { klasse: sfk },
        schraubengroesse: 36.0,
        gewindesteigung: 4.,
        passschraube,
      },
      _ => panic!(),
    }
  }
  pub fn get_schrauben_schneider() -> [u32; 8] {
    [12, 16, 20, 22, 24, 27, 30, 36]
  }

  // DIN 13-1 ISO Schrauben!
  fn flankendurchmesser(&self) -> f64 {
    let h = 3.0_f64.sqrt() / 2.0 * self.gewindesteigung;
    let dp = self.schraubengroesse - 2.0 * 3. / 8. * h;
    dp
  }
  // in mm
  fn kerndurchmesser(&self) -> f64 {
    let h = 3.0_f64.sqrt() / 2.0 * self.gewindesteigung;
    let D2 = self.schraubengroesse - 3. / 4. * h;
    D2 - 2. * (h / 2. - h / 6.)
  }
  /// in mm
  fn nennlochdurchmesser(&self) -> f64 {
    match self.schraubengroesse as u32 {
      12 => self.schraubengroesse + 1.0,
      16 => self.schraubengroesse + 2.0,
      20 => self.schraubengroesse + 2.0,
      22 => self.schraubengroesse + 2.0,
      24 => self.schraubengroesse + 2.0,
      27 => self.schraubengroesse + 3.0,
      30 => self.schraubengroesse + 3.0,
      36 => self.schraubengroesse + 3.0,
      _ => panic!(),
    }
  }

  /// in `mm`
  fn schluesselweite(&self) -> f64 {
    match self.schraubengroesse as u32 {
      12 => 18.,
      16 => 24.,
      20 => 30.,
      22 => 34.,
      24 => 36.,
      27 => 41.,
      30 => 46.,
      36 => 55.,
      _ => panic!(),
    }
  }

  /// in `mm`
  fn eckenmass(&self) -> f64 {
    self.schluesselweite() * 2.0 / (3.0_f64).sqrt()
  }

  /// in `mm`
  pub fn mittlerer_schraubenkopfdurchmesser(&self) -> f64 {
    return self.schluesselweite() * 0.5 + self.eckenmass() * 0.5;
  }

  /// in `mm^2`
  pub fn spannungsflaeche_schraube(&self) -> f64 {
    ((self.flankendurchmesser() + self.kerndurchmesser()) / 2.).powi(2) / 4. * consts::PI
  }

  pub fn schaftdurchmesser(&self) -> f64 {
    if self.passschraube {
      self.schraubengroesse + 1.0
    } else {
      self.schraubengroesse
    }
  }

  /// in `mm^2`
  pub fn flaeche_schraube(&self) -> f64 {
    let d = self.schaftdurchmesser();
    d.powi(2) / 4.0 * consts::PI
  }

  /// in `N`
  pub fn abscherkraft_schraube(&self, ist_gewinde: bool) -> f64 {
    if !ist_gewinde {
      0.6 * self.festigkeitsklasse.zugfestigkeit() * self.flaeche_schraube() / 1.25
    } else if self.festigkeitsklasse.klasse == [19, 9] {
      0.5 * self.festigkeitsklasse.zugfestigkeit() * self.spannungsflaeche_schraube() / 1.25
    } else {
      0.6 * self.festigkeitsklasse.zugfestigkeit() * self.spannungsflaeche_schraube() / 1.25
    }
  }

  /// in `N`
  pub fn zugkraft_schraube(&self, ist_senkschraube: bool) -> f64 {
    let k = if ist_senkschraube { 0.63 } else { 0.9 };
    return k * self.festigkeitsklasse.zugfestigkeit() * self.spannungsflaeche_schraube() / 1.25;
  }
}

impl ISOSchraube {
  pub fn durchstanzwiderstand(&self, blechdicke: f64, zugfestigkeit_blech: f64) -> f64 {
    return 0.6
      * consts::PI
      * self.mittlerer_schraubenkopfdurchmesser()
      * blechdicke
      * zugfestigkeit_blech
      / 1.25;
  }

  pub fn lochleibungswiderstand(
    &self,
    blechdicke: f64,
    zugfestigkeit_blech: f64,
    p1: f64,
    e1: f64,
    p2: f64,
    e2: f64,
    innenliegend: bool,
  ) -> f64 {
    let alpha1;
    let k_1;
    let d0 = self.nennlochdurchmesser();
    if innenliegend {
      alpha1 = (p1 / (3.0 * d0) - 0.25)
        .min(self.festigkeitsklasse.zugfestigkeit() / zugfestigkeit_blech)
        .min(1.);
      k_1 = (1.4 * p2 / d0 - 1.7).min(2.5);
    } else {
      alpha1 = (e1 / (3. * d0))
        .min(self.festigkeitsklasse.zugfestigkeit() / zugfestigkeit_blech)
        .min(1.);
      k_1 = (2.8 * e2 / (d0) - 1.7).min(1.4 * p2 / d0 - 1.7).min(2.5)
    }

    return k_1 * alpha1 * zugfestigkeit_blech * self.schaftdurchmesser() * blechdicke / 1.25;
  }
}
