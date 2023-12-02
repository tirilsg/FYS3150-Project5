#include "diffraction.hpp"
#include <armadillo>
#include <iostream>
#include <complex>
#include <fstream>
#include <assert.h>


Diff::Diff(double h_in, double dt_in, double T_in){
    h = h_in;
    dt = dt_in;
    T = T_in;
    M = 1. / h + 1;
    N = T / dt + 1;
    r = arma::cx_double(0, dt / (2 * h * h));
    U.zeros(M - 2, M - 2);
    V.zeros(M - 2, M - 2);
    S.zeros(M - 2, M - 2, N);
    A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
}

int Diff::index_k(int i, int j){
  return i + (M - 2) * j;
}

void Diff::set_potential(int slits, int V0){
  V.zeros(M - 2, M - 2);
  double length = 0.05;
  double aperture = 0.05;
  double thick = 0.02;
  int c_x = (M - 2) / 2; 
  int c_y = (M - 2) / 2; 
  int c_l = static_cast<int>(length / h); 
  int a = static_cast<int>(aperture / h);     
  int t = static_cast<int>(thick / h);   
  int n_c = slits - 1; 
  int e_l = (M - 2) - slits * a - n_c * c_l;
  e_l /= 2;
  for (int i = 0; i < (M - 2); i++){
    for (int j = c_y - t / 2; j <= c_y + t / 2; j++){
      if (slits == 1){
        if (i < e_l || i > e_l + a){
          V(i, j) = V0;
        }
      }
      else if (slits == 2){
        if (i < e_l || (i >= e_l + a && i <= e_l + a + c_l) || i > e_l + 2 * a + c_l){
          V(i, j) = V0;
        }
      }
      else if (slits == 3){
        if (i < e_l || (i >= e_l + a && i <= e_l + a + c_l) ||(i > e_l + 2 * a + c_l && i <= e_l + 2 * (a + c_l)) || i > e_l + 2 * a + 3 * c_l){
          V(i, j) = V0;
        }
      }
    }
  }
}

void Diff::fill_matrices(){
  arma::cx_vec a((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_vec b((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_mat sub_diag(M - 2, M - 2, arma::fill::zeros);
  sub_diag.diag(-1).fill(r);
  sub_diag.diag(1).fill(r);
  A.diag(M - 2).fill(-r);
  A.diag(2 - M).fill(-r);
  B.diag(M - 2).fill(r);
  B.diag(2 - M).fill(r);
  arma::cx_double a_k, b_k;
  int k;
  for (int i = 0; i < M - 2; i++){
    for (int j = 0; j < M - 2; j++){
      k = index_k(i, j);
      a_k = arma::cx_double(1., 4. * r.imag() + .5 * dt * V(i, j));
      b_k = std::conj(a_k);
      a(k) = a_k;
      b(k) = b_k;
    }
  }
  for (int i = 0; i < M - 2; i++){
    int j = i * (M - 2);
    A.submat(j, j, j + M - 3, j + M - 3) = -sub_diag;
    B.submat(j, j, j + M - 3, j + M - 3) = sub_diag;
  }
  A.diag() = a;
  B.diag() = b;
}

void Diff::simulate(){
  for (int i = 1; i < N; i++){
    u_current = U.as_col();
    u_new = arma::spsolve(A, B * u_current);
    for (int n = 0; n < U.n_cols; n ++){
      for (int k = 0; k < M - 2; k++){
        U.col(n)(k) = u_new(k + n * (M - 2));
      }
    }
    S.slice(i) = U;
    u_current = u_new;
  }
}

void Diff::set_initial_state(double xc, double sigmax, double px, double yc, double sigmay, double py){
  double re, im;
  int xi, yi;
  for (int i = 0; i < M - 2; i++){
    yi = i;
    for (int j = 0; j < M - 2; j++){
      xi = j;
      re = -(xi * h - xc) * (xi * h * xc) / (2 * sigmax * sigmax)-(yi * h - yc) * (yi * h - yc) / (2 * sigmay * sigmay);
      im = px * (xi * h - xc) + (yi * h - yc);
      U(i, j) = std::exp(arma::cx_double(re, im));
    }
  }
  U /= std::sqrt(arma::accu(arma::conj(U) % U));
  S.slice(0) = U;
}










/*Diff::Diff(double h_in, double dt_in, double T_in){
    h = h_in;
    dt = dt_in;
    T = T_in;
    M = 1. / h + 1;
    N = T / dt + 1;
    r = arma::cx_double(0, dt / (2 * h * h));
    U.zeros(M - 2, M - 2);
    V.zeros(M - 2, M - 2);
    S.zeros(M - 2, M - 2, N);
    A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
}

int Diff::index_k(int i, int j){
  return i + (M - 2) * j;
}

void Diff::set_potential(int slits, int V0){
  V.zeros(M - 2, M - 2);
  double length = 0.05;
  double aperture = 0.05;
  double thick = 0.02;
  int c_x = (M - 2) / 2; 
  int c_y = (M - 2) / 2; 
  int c_l = static_cast<int>(length / h); 
  int a = static_cast<int>(aperture / h);     
  int t = static_cast<int>(thick / h);   
  int n_c = slits - 1; 
  int e_l = (M - 2) - slits * a - n_c * c_l;
  e_l /= 2;
  for (int i = 0; i < (M - 2); i++){
    for (int j = c_y - t / 2; j <= c_y + t / 2; j++){
      if (slits == 1){
        if (i < e_l || i > e_l + a){
          V(i, j) = V0;
        }
      }
      else if (slits == 2){
        if (i < e_l || (i >= e_l + a && i <= e_l + a + c_l) || i > e_l + 2 * a + c_l){
          V(i, j) = V0;
        }
      }
      else if (slits == 3){
        if (i < e_l || (i >= e_l + a && i <= e_l + a + c_l) ||(i > e_l + 2 * a + c_l && i <= e_l + 2 * (a + c_l)) || i > e_l + 2 * a + 3 * c_l){
          V(i, j) = V0;
        }
      }
    }
  }
}

void Diff::fill_matrices(){
  arma::cx_vec a((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_vec b((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_mat sub_diag(M - 2, M - 2, arma::fill::zeros);
  sub_diag.diag(-1).fill(r);
  sub_diag.diag(1).fill(r);
  A.diag(M - 2).fill(-r);
  A.diag(2 - M).fill(-r);
  B.diag(M - 2).fill(r);
  B.diag(2 - M).fill(r);
  arma::cx_double a_k, b_k;
  int k;
  for (int i = 0; i < M - 2; i++){
    for (int j = 0; j < M - 2; j++){
      k = index_k(i, j);
      a_k = arma::cx_double(1., 4. * r.imag() + .5 * dt * V(i, j));
      b_k = std::conj(a_k);
      a(k) = a_k;
      b(k) = b_k;
    }
  }
  for (int i = 0; i < M - 2; i++){
    int j = i * (M - 2);
    A.submat(j, j, j + M - 3, j + M - 3) = -sub_diag;
    B.submat(j, j, j + M - 3, j + M - 3) = sub_diag;
  }
  A.diag() = a;
  B.diag() = b;
}

void Diff::simulate(){
  for (int i = 1; i < N; i++){
    u_current = U.as_col();
    u_new = arma::spsolve(A, B * u_current);
    for (int n = 0; n < U.n_cols; n ++){
      for (int k = 0; k < M - 2; k++){
        U.col(n)(k) = u_new(k + n * (M - 2));
      }
    }
    S.slice(i) = U;
    u_current = u_new;
  }
}

void Diff::set_initial_state(double xc, double sigmax, double px, double yc, double sigmay, double py){
  //arma::cx_double exponent;
  double re, im;
  int xi, yi;
  for (int i = 0; i < M - 2; i++){
    yi = i;
    for (int j = 0; j < M - 2; j++){
      xi = j;
      re = -(xi * h - xc) * (xi * h * xc) / (2 * sigmax * sigmax)-(yi * h - yc) * (yi * h - yc) / (2 * sigmay * sigmay);
      im = px * (xi * h - xc) + (yi * h - yc);
      U(i, j) = std::exp(arma::cx_double(re, im));
    }
  }
  U /= std::sqrt(arma::accu(arma::conj(U) % U));
  S.slice(0) = U;
}*/









///////////////////////////////////////////////NOTES

/*Diff::Diff(double h_in, double dt_in, double T_in){
    h = h_in;
    dt = dt_in;
    T = T_in;
    M = 1. / h + 1;
    N = T / dt + 1;
    r = arma::cx_double(0, dt / (2 * h * h));
    U.zeros(M - 2, M - 2);
    V.zeros(M - 2, M - 2);
    S.zeros(M - 2, M - 2, N);
    A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
}

int Diff::index_k(int i, int j){
  return i + (M - 2) * j;
}

void Diff::set_potential(int slits){
  V.zeros(M - 2, M - 2);
  double center_length = 0.05;
  double aperture = 0.05;
  double thickness = 0.02;
  int c_x = (M - 2) / 2; 
  int c_y = (M - 2) / 2; 
  int c_l = static_cast<int>(center_length / h); 
  int a = static_cast<int>(aperture / h);     
  int t = static_cast<int>(thickness / h);   
  int n_c = slits - 1; 
  int e_l = (M - 2) - slits * a - n_c * c_l;
  e_l /= 2;
  for (int i = 0; i < (M - 2); i++){
    for (int j = c_y - t / 2; j <= c_y + t / 2; j++){
      if (slits == 1){
        if (i < e_l || i > e_l + a){
          V(i, j) = 1e10;
        }
      }
      else if (slits == 2){
        if (i < e_l || (i >= e_l + a && i <= e_l + a + c_l) || i > e_l + 2 * a + c_l){
          V(i, j) = 1e10;
        }
      }
      else if (slits == 3){
        if (i < e_l || (i >= e_l + a && i <= e_l + a + c_l) ||(i > e_l + 2 * a + c_l && i <= e_l + 2 * (a + c_l)) || i > e_l + 2 * a + 3 * c_l){
          V(i, j) = 1e10;
        }
      }
    }
  }
}

void Diff::fill_matrices(){
  arma::cx_vec a((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_vec b((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_mat sub_diag(M - 2, M - 2, arma::fill::zeros);
  sub_diag.diag(-1).fill(r);
  sub_diag.diag(1).fill(r);
  A.diag(M - 2).fill(-r);
  A.diag(2 - M).fill(-r);
  B.diag(M - 2).fill(r);
  B.diag(2 - M).fill(r);
  arma::cx_double a_k, b_k;
  int k;
  for (int i = 0; i < M - 2; i++){
    for (int j = 0; j < M - 2; j++){
      k = index_k(i, j);
      a_k = arma::cx_double(1., 4. * r.imag() + .5 * dt * V(i, j));
      b_k = std::conj(a_k);
      a(k) = a_k;
      b(k) = b_k;
    }
  }
  for (int i = 0; i < M - 2; i++){
    int j = i * (M - 2);
    A.submat(j, j, j + M - 3, j + M - 3) = -sub_diag;
    B.submat(j, j, j + M - 3, j + M - 3) = sub_diag;
  }
  A.diag() = a;
  B.diag() = b;
}

void Diff::solve(){
  for (int i = 1; i < N; i++){
    u_current = U.as_col();
    u_new = arma::spsolve(A, B * u_current);
    for (int n = 0; n < U.n_cols; n ++){
      for (int k = 0; k < M - 2; k++){
        U.col(n)(k) = u_new(k + n * (M - 2));
      }
    }
    S.slice(i) = U;
    u_current = u_new;
  }
}

void Diff::set_initial_state(double x_c, double sigma_x, double p_x, double y_c, double sigma_y, double p_y){
  arma::cx_double exponent;
  double re, im;
  int xi, yi;
  for (int i = 0; i < M - 2; i++){
    yi = i;
    for (int j = 0; j < M - 2; j++){
      xi = j;
      re = -(xi * h - x_c) * (xi * h * x_c) / (2 * sigma_x * sigma_x)-(yi * h - y_c) * (yi * h - y_c) / (2 * sigma_y * sigma_y);
      im = p_x * (xi * h - x_c) + (yi * h - y_c);
      exponent = arma::cx_double(re, im);
      U(i, j) = std::exp(exponent);
    }
  }
  U /= std::sqrt(arma::accu(arma::conj(U) % U));
  S.slice(0) = U;
}*/



/*Diff::Diff(double r_, double dt_, double L_, double W_, double V0_, double b_, double h_, double m_, double x0_, double sigma_, int M_, double x_max_, double y_max_, double t_max_, double sigma_t_,
           double n_slit_, double thickx_, double centerx_, double slit_sep_, double aperture_, double xc_, double yc_, double widthx_, double widthy_, double px_, double py_)
    : r(r_), dt(dt_), L(L_), W(W_), V0(V0_), b(b_), h(h_), m(m_), x0(x0_), sigma(sigma_), M(M_), x_max(x_max_), y_max(y_max_), t_max(t_max_), sigma_t(sigma_t_),
      n_slit(n_slit_), thickx(thickx_), centerx(centerx_), slit_sep(slit_sep_), aperture(aperture_), xc(xc_), yc(yc_), widthx(widthx_), widthy(widthy_), px(px_), py(py_) {
    // Initialize other members and perform any necessary setup
    M = 1. / h + 1;
    len = (M - 2);
    V = potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    u = wave_init(xc, yc, widthx, widthy, px, py, M);
    create_AB(V, h, dt, M);
    A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    r = arma::cx_double(0, dt / (2 * h * h));
    n_timesteps = std::round(t_max / dt);
    storage = arma::cube(len, len, n_timesteps + 1);
    Re = arma::cube(len, len, n_timesteps + 1);
    Im = arma::cube(len, len, n_timesteps + 1);
    T = t_max;
}


arma::mat Diff::potential(double potential, int M, double n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    arma::mat V = arma::zeros<arma::mat>(M, M);

    // Simple example: Potential well
    if (potential == 1) {
        // Set the potential well inside a specified region
        int start_x = static_cast<int>((centerx - thickx / 2.0) / h);
        int end_x = static_cast<int>((centerx + thickx / 2.0) / h);
        int start_y = static_cast<int>((n_slit * slit_sep - aperture / 2.0) / h);
        int end_y = static_cast<int>((n_slit * slit_sep + aperture / 2.0) / h);

        V.submat(start_x, start_y, end_x, end_y).fill(V0);
    }

    // Add more cases for different potentials if needed

    return V;
}
arma::cx_mat Diff::wave_init(double xc, double yc, double widthx, double widthy, double px, double py, int M) {
    arma::cx_mat u = arma::zeros<arma::cx_mat>(M, M);

    // Simple example: Gaussian wave packet
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            double x = i * h;
            double y = j * h;

            double exponent = -0.5 * ((x - xc) * (x - xc) / (widthx * widthx) + (y - yc) * (y - yc) / (widthy * widthy));
            u(i, j) = std::exp(arma::cx_double(0, (px * (x - xc) + py * (y - yc)))) * std::exp(exponent);
        }
    }

    return u;
}*/

//////////////////////herfra
/*int Diff::index_k(int i, int j)
{
  return i + (M - 2) * j;
}


int Diff::pair_to_single(const int i, const int j, const int len) {
    return (j) * len + (i);
}

std::tuple<int, int> Diff::single_to_pair(const int k, const int len) {
    int i = k % len;
    int j = k / len;
    return std::tuple<int, int>{i, j};
}

std::vector<arma::sp_cx_dmat> Diff::create_AB(arma::mat &V, const double h, const double dt, const int M) {
    int len = M-2;
    std::complex<double> im(0., 1.);
    std::complex<double> r = im*dt/(2*h*h);
    arma::cx_dvec a = arma::cx_dvec(len*len);
    arma::cx_dvec b = arma::cx_dvec(len*len);
    for (int ii=0; ii < (len); ii++) {
        arma::vec temp_col = V.col(ii);
        for (int i=0; i < (len); i++) {
            std::complex<double> temp_im = im * ((dt*temp_col(i)) / 2.);
            a(pair_to_single(i, ii, len)) = 1. + 4.*r + temp_im;
            b(pair_to_single(i, ii, len)) = 1. - 4.*r - temp_im;
        }
    }
    return std::vector<arma::sp_cx_dmat>{create_mat(a, -r, len), create_mat(b, r, len)};
}

arma::sp_cx_dmat Diff::create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len) {
    int lenlen = len*len;
    arma::sp_cx_dmat A = arma::sp_cx_dmat(lenlen, lenlen);
    for (int i = 0; i < len; i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start, start, end, end) = create_tri(a, r, len, i);
    }
    for (int i = 1; i<(len); i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start-len, start, end-len, end) = create_rdiag(r, len);
        A.submat(start, start-len, end, end-len) = create_rdiag(r, len);
    }
    return A;
}

arma::sp_cx_dmat Diff::create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    for (int ii = 0; ii<len; ii++) {
        temp(ii, ii) = a(pair_to_single(ii, i, len));
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }
    return temp;
}

arma::sp_cx_dmat Diff::create_rdiag(const std::complex<double> r, const int len) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    temp.diag().fill(r);
    return temp;
}*/


/*Diff::Diff(double r_, double dt_, double L_, double W_, double V0_, double b_, double h_, double m_, double x0_, double sigma_, int M_, double x_max_, double y_max_, double t_max_, double sigma_t_)
    : r(r_), dt(dt_), L(L_), W(W_), V0(V0_), b(b_), h(h_), m(m_), x0(x0_), sigma(sigma_), M(M_), x_max(x_max_), y_max(y_max_), t_max(t_max_), sigma_t(sigma_t_) {
    create_AB(V, h, dt, M);
    len = (M - 2);
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    this->u = wave_init(xc, yc, widthx, widthy, px, py, M);
    // Assuming create_AB modifies AB directly, or you get the necessary information from create_AB
    A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    r = arma::cx_double(0, dt / (2 * h * h));
    n_timesteps = std::round(T/dt);
    storage = arma::cube(len, len, n_timesteps + 1);
    Re = arma::cube(len, len, n_timesteps + 1);
    Im = arma::cube(len, len, n_timesteps + 1);
    this->T = T;
    this->dt = dt;
}*/


/*Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    M = std::round(1./h);
    len = (M-2);
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    this->u = wave_init(xc, yc, widthx, widthy, px, py, M);
    AB = create_AB(V, h, dt, M);
    n_timesteps = std::round(T/dt);
    storage = arma::cube(len, len, n_timesteps + 1);
    Re = arma::cube(len, len, n_timesteps + 1);
    Im = arma::cube(len, len, n_timesteps + 1);
    this->T = T;
    this->dt = dt;
}

void Diff::run() {
    storage.slice(0) = this->probability(this->u);
    print_u(0);
    for (int i=1; i < n_timesteps + 1 ; i++) {
        b = AB.at(1) * this->u;
        this->u = arma::spsolve(AB.at(0), b);
        storage.slice(i) = this->probability(this->u);
        print_u(i);
    }
}

arma::mat Diff::probability(arma::cx_dvec &u) {
    arma::mat prob = arma::mat(len, len);
    for (int j=0; j < len; j++) {
        for (int i=0; i < len; i++) {
            prob(i, j) = std::real(std::conj(this->u(Diff::pair_to_single(i, j, len))) * this->u(Diff::pair_to_single(i, j, len)));
        }
    }
    return prob;
}

void Diff::print_u(int t) {

    arma::cx_mat temp = arma::cx_mat(len,len);
    int temp_len = len*len;
    for (int i = 0; i < temp_len; i++) {
        std::tuple<int,int> a = single_to_pair(i,len);
        temp(std::get<0>(a),std::get<1>(a)) = this->u(i);
        }
    Re.slice(t) = arma::real(temp);
    Im.slice(t) = arma::imag(temp);

}

arma::cx_dvec Diff::wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
    int len = M-2;
    arma::cx_dvec u = arma::cx_dvec(len*len);
    double h = 1./len;
    std::complex<double> im(0., 1.);
    for (int i=0; i<(len); i++) {
        for (int j=0; j<(len); j++) {
            double temp_x = ((h*(j+1)) - centerx);
            double temp_y = ((h*(i+1)) - centery);
            std::complex<double> arg = std::exp(-(temp_x*temp_x/(2*widthx*widthx)) - (temp_y*temp_y/(2*widthy*widthy)) + im*(px*temp_x) + im*(py*temp_y));
            u(pair_to_single(i, j, len)) = arg;
        }
    }
    u = u / std::sqrt(arma::accu(arma::conj(u) % u));

    return u;
}



arma::mat Diff::potential(double potential, int M, int slits, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M - 2;
    double h = 1.0 / len;
    arma::mat V = arma::zeros<arma::mat>(len, len);

    int c_x = len / 2; // Center index along x axis
    int c_y = len / 2; // Center index along y axis
    int c_l = std::round(centerx * len); // Index length of center wall
    int a = std::round(aperture * len); // Index length of aperture
    int t = std::round(thickx * len);   // Index length of wall thickness

    int n_c = slits - 1; // Number of center pieces
    int e_l = len - slits * a - n_c * c_l; // Index length of end wall pieces
    e_l /= 2;

    for (int i = 0; i < len; ++i) {
        for (int j = c_y - t / 2; j <= c_y + t / 2; ++j) {

            if (slits == 1) {
                if (i < e_l || e_l + a < i) {
                    V(i, j) = potential;
                }
            }

            if (slits == 2) {
                if (i < e_l || (e_l + a <= i && i <= e_l + a + c_l) || e_l + 2 * a + c_l < i) {
                    V(i, j) = potential;
                }
            }

            if (slits == 3) {
                if (i < e_l || (e_l + a <= i && i <= e_l + a + c_l) || (e_l + 2 * a + c_l < i && i <= e_l + 2 * (a + c_l)) || e_l + 2 * a + 3 * c_l < i) {
                    V(i, j) = potential;
                }
            }
        }
    }

    return V;
}*/

/*arma::mat Diff::potential(double potential, int M, int slits, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M - 2;
    double h = 1.0 / len;
    arma::mat V = arma::mat(len, len, arma::fill::zeros);
    int c_x = len / 2;                               // Center index along x axis
    int c_y = len / 2;                               // Center index along y axis
    int c_l = static_cast<int>(centerx * len);       // Index length of center wall
    int a = static_cast<int>(aperture * len);        // Index length of aperture
    int t = static_cast<int>(thickx * len);          // Index length of wall thickness
    int n_c = slits - 1;                           // Number of center pieces
    int e_l = len - slits * a - n_c * c_l;         // Index length of end wall pieces
    e_l /= 2;
    for (int i = 0; i < len; ++i) {
        for (int j = c_y - t / 2; j <= c_y + t / 2; ++j) {
            if (slits == 1) {
                if (i < e_l || e_l + a < i) {
                    V(i, j) = potential;
                }
            } else if (slits == 2) {
                if (i < e_l || (e_l + a <= i && i <= e_l + a + c_l) || e_l + 2 * a + c_l < i) {
                    V(i, j) = potential;
                }
            } else if (slits == 3) {
                if (i < e_l || (e_l + a <= i && i <= e_l + a + c_l) || (e_l + 2 * a + c_l < i && i <= e_l + 2 * (a + c_l)) || e_l + 2 * a + 3 * c_l < i) {
                    V(i, j) = potential;
                }
            }
        }
    }
    return V;
}*/
/*arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M - 2;
    double h = 1.0 / len;
    arma::mat V = arma::mat(len, len, arma::fill::zeros);
    if (n_slit > 0) {
        int slitWidth = std::round(thickx / h);
        int slitSeparation = std::round(slit_sep / h);
        int start = std::round((centerx - 0.5 * thickx) / h);
        for (int i = 0; i < n_slit; ++i) {
            for (int j = start; j < start + slitWidth; ++j) {
                V.col(j).fill(potential);
            }
            start += slitWidth + slitSeparation;
        }
    }
    return V;
}*/

/*arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M-2;
    double h = 1./len;
    arma::mat V = arma::mat(len, len).zeros();
    int j_start = std::round((centerx - (thickx/2))/h);
    int j_end = std::round((centerx + (thickx/2))/h) - 1;
    bool center_open = n_slit%2 != 0;   
    int counter = n_slit;
    int i_center = std::round(0.5/h) -1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture/h);
    int len_slit_sep = std::round(slit_sep/h);
    if (center_open) {
        int i_loc_up = i_center - std::round(aperture/(2*h));
        int i_loc_down = i_center + std::round(aperture/(2*h));
        i_loc_up--;
        i_loc_down++;
        counter--; 
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }
    else if (n_slit == 0) {
        for (int j=j_start; j<j_end; j++) {
            V.col(j).fill(potential);
        }
        return V;
    }
    else {
        int temp = std::round(len_slit_sep/2.);
        for (int ii=0; ii<temp; ii++) {
            for (int j=j_start; j<j_end; j++) {
                V(i_loc_up, j) = potential;
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }
}*/




/*arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M-2;
    double h = 1./len;
    arma::mat V = arma::mat(len, len).zeros();
    int j_start = std::round((centerx - (thickx/2))/h);
    int j_end = std::round((centerx + (thickx/2))/h) - 1;
    bool center_open = n_slit%2 != 0;   
    int counter = n_slit;
    int i_center = std::round(0.5/h) -1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture/h);
    int len_slit_sep = std::round(slit_sep/h);
    if (center_open) {
        int i_loc_up = i_center - std::round(aperture/(2*h));
        int i_loc_down = i_center + std::round(aperture/(2*h));
        i_loc_up--;
        i_loc_down++;
        counter--; 
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }
    else if (n_slit == 0) {
        for (int j=j_start; j<j_end; j++) {
            V.col(j).fill(potential);
        }
        return V;
    }
    else {
        int temp = std::round(len_slit_sep/2.);
        for (int ii=0; ii<temp; ii++) {
            for (int j=j_start; j<j_end; j++) {
                V(i_loc_up, j) = potential;
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }
}*/



/*Matrix::Matrix(double h_in, double dt_in, double T_in, double x_c, double sigma_x, double p_x, double y_c,
               double sigma_y, double p_y, int slits)
{
    h = h_in;
    dt = dt_in;
    T = T_in;

    M = 1. / h + 1;
    N = T / dt + 1;
    r = arma::cx_double(0, dt / (2 * h * h));

    U.zeros(M - 2, M - 2);
    V.zeros(M - 2, M - 2);
    S.zeros(M - 2, M - 2, N);

    A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));

    set_potential(slits);
    fill_matrices();
    set_initial_state(x_c, sigma_x, p_x, y_c, sigma_y, p_y);
}

int Matrix::index_k(int i, int j)
{
    return i + (M - 2) * j;
}

void Matrix::set_potential(int slits)
{
    V.zeros();

    double V_0 = 1e10;          // Potential amplitude, adjust as needed
    double thickness = 0.02;    // Adjust the thickness of the walls
    double center_length = 0.05; // Adjust the length of the center wall
    double apperture = 0.05;    // Adjust the length of the apperture

    int c_x = (M - 2) / 2;                          // Center index along x axis
    int c_y = (M - 2) / 2;                          // Center index along y axis
    int c_l = static_cast<int>(center_length / h); // Index length of center wall
    int a = static_cast<int>(apperture / h);       // Index length of apperture
    int t = static_cast<int>(thickness / h);

    int n_c = slits - 1; // Number of center pieces
    int e_l = (M - 2) - slits * a - n_c * c_l; // Index length of end wall pieces
    e_l /= 2;

    for (int i = 0; i < M - 2; ++i)
    {
        for (int j = c_y - t / 2; j <= c_y + t / 2; ++j)
        {
            if (slits == 0)
            {
                // No slits, potential is zero
                V(i, j) = 0.0;
            }
            else if (slits == 1)
            {
                if (i < e_l || e_l + a < i)
                    V(i, j) = V_0;
            }
            else if (slits == 2)
            {
                if (i < e_l || (e_l + a <= i && i <= e_l + a + c_l) || e_l + 2 * a + c_l < i)
                    V(i, j) = V_0;
            }
            else if (slits == 3)
            {
                if (i < e_l || (e_l + a <= i && i <= e_l + a + c_l) || (e_l + 2 * a + c_l < i && i <= e_l + 2 * (a + c_l)) || e_l + 2 * a + 3 * c_l < i)
                    V(i, j) = V_0;
            }
            // Add more cases as needed for different slit configurations
        }
    }
}

void Matrix::fill_matrices()
{
    A.zeros();
    B.zeros();
    arma::cx_vec a((M - 2) * (M - 2), arma::fill::zeros);
    arma::cx_vec b((M - 2) * (M - 2), arma::fill::zeros);
    arma::cx_mat sub_diag(M - 2, M - 2, arma::fill::zeros);
    sub_diag.diag(-1).fill(r);
    sub_diag.diag(1).fill(r);
    A.diag(M - 2).fill(-r);
    A.diag(2 - M).fill(-r);
    B.diag(M - 2).fill(r);
    B.diag(2 - M).fill(r);
    arma::cx_double a_k, b_k;
    int k;
    for (int i = 0; i < M - 2; i++){
        for (int j = 0; j < M - 2; j++){
            k = index_k(i, j);
            a_k = arma::cx_double(1., 4. * r.imag() + .5 * dt * V(i, j));
            b_k = std::conj(a_k);
            a(k) = a_k;
            b(k) = b_k;
        }
    }
    for (int i = 0; i < M - 2; i++){
        int j = i * (M - 2);
        A.submat(j, j, j + M - 3, j + M - 3) = -sub_diag;
        B.submat(j, j, j + M - 3, j + M - 3) = sub_diag;
    }
    A.diag() = a;
    B.diag() = b;
}

void Matrix::solve(){
    for (int i = 1; i < N; i++){
        u_current = U.as_col();
        u_new = arma::spsolve(A, B * u_current);
        for (int n = 0; n < U.n_cols; n++){
            for (int k = 0; k < M - 2; k++){
                U.col(n)(k) = u_new(k + n * (M - 2));
            }
        }
        S.slice(i) = U;
        u_current = u_new;
    }
}

void Matrix::set_initial_state(double x_c, double sigma_x, double p_x, double y_c, double sigma_y, double p_y){
    arma::cx_double exponent;
    double re, im;
    int xi, yi;
    for (int i = 0; i < M - 2; i++){
        yi = i;
        for (int j = 0; j < M - 2; j++){
            xi = j;
            re = -(xi * h - x_c) * (xi * h - x_c) / (2 * sigma_x * sigma_x)-(yi * h - y_c) * (yi * h - y_c) / (2 * sigma_y * sigma_y);
            im = p_x * (xi * h - x_c) + p_y * (yi * h - y_c);
            exponent = arma::cx_double(re, im);
            U(i, j) = std::exp(exponent);
        }
    }
    U /= std::sqrt(arma::accu(arma::conj(U) % U));
    S.slice(0) = U;
}*/



/*Matrix::Matrix(double h_in, double dt_in, double T_in)
{
  h = h_in;
  dt = dt_in;
  T = T_in;

  M = 1. / h + 1;
  N = T / dt + 1;
  r = arma::cx_double(0, dt / (2 * h * h));

  U.zeros(M - 2, M - 2);
  V.zeros(M - 2, M - 2);
  S.zeros(M - 2, M - 2, N);

  A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
  B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
}

// Get index k from (i, j)
int Matrix::index_k(int i, int j)
{
  return i + (M - 2) * j;
}

// Set potential barriers
void Matrix::set_potential(std::string potential)
{
  // Load the file into the potential
  V.load(potential, arma::raw_ascii);
  assert(V.n_rows == (M - 2) && V.n_cols == (M - 2));

}

// Fill the matrices according to the Crank-Nicolson regime
void Matrix::fill_matrices()
{
  // Vectors a and b, which will be the main diagonals of A and B
  arma::cx_vec a((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_vec b((M - 2) * (M - 2), arma::fill::zeros);

  // Tridiagonal sub matrix with signature (r, 0, r)
  arma::cx_mat sub_diag(M - 2, M - 2, arma::fill::zeros);
  sub_diag.diag(-1).fill(r);
  sub_diag.diag(1).fill(r);

  // Filling diagonals of A and B
  A.diag(M - 2).fill(-r);
  A.diag(2 - M).fill(-r);
  B.diag(M - 2).fill(r);
  B.diag(2 - M).fill(r);

  // Filling a and b
  arma::cx_double a_k, b_k;
  int k;

  for (int i = 0; i < M - 2; i++)
  {
    for (int j = 0; j < M - 2; j++)
    {
      k = index_k(i, j);

      a_k = arma::cx_double(1., 4. * r.imag() + .5 * dt * V(i, j));
      b_k = std::conj(a_k);

      a(k) = a_k;
      b(k) = b_k;
    }
  }

  // Filling A and B with sub matrices
  for (int i = 0; i < M - 2; i++)
  {
    int j = i * (M - 2);

    A.submat(j, j, j + M - 3, j + M - 3) = -sub_diag;
    B.submat(j, j, j + M - 3, j + M - 3) = sub_diag;
  }

  // Setting main diagonals of A and B
  A.diag() = a;
  B.diag() = b;
}

// Solve matrix eq. Au^(n+1) = Bu^n
void Matrix::solve()
{
  // Solving the equation for each time step and saving as slice in S
  for (int i = 1; i < N; i++)
  {
    u_current = U.as_col();
    u_new = arma::spsolve(A, B * u_current);

    // Refilling u vector into U matrix
    for (int n = 0; n < U.n_cols; n ++)
    {
      for (int k = 0; k < M - 2; k++)
      {
        U.col(n)(k) = u_new(k + n * (M - 2));
      }
    }

    S.slice(i) = U;
    u_current = u_new;
  }
}

// Set the initial state of the system
void Matrix::set_initial_state(double x_c, double sigma_x, double p_x, double y_c,
                       double sigma_y, double p_y)
{
  arma::cx_double exponent;
  double re, im;
  int xi, yi;

  for (int i = 0; i < M - 2; i++)
  {
    yi = i;

    for (int j = 0; j < M - 2; j++)
    {
      xi = j;

      re = -(xi * h - x_c) * (xi * h * x_c) / (2 * sigma_x * sigma_x)
           -(yi * h - y_c) * (yi * h - y_c) / (2 * sigma_y * sigma_y);

      im = p_x * (xi * h - x_c) + (yi * h - y_c);

      exponent = arma::cx_double(re, im);

      U(i, j) = std::exp(exponent);
    }
  }

  // Normalize
  U /= std::sqrt(arma::accu(arma::conj(U) % U));
  S.slice(0) = U;
}*/












/*int Diff::pair_to_single(const int i, const int j, const int len) {
    return (j) * len + (i);
}

std::tuple<int, int> Diff::single_to_pair(const int k, const int len) {
    int i = k % len;
    int j = k / len;
    return std::tuple<int, int>{i, j};
}

std::vector<arma::sp_cx_dmat> Diff::create_AB(arma::mat &V, const double h, const double dt, const int M) {
    int len = M-2;
    std::complex<double> im(0., 1.);
    std::complex<double> r = im*dt/(2*h*h);
    arma::cx_dvec a = arma::cx_dvec(len*len);
    arma::cx_dvec b = arma::cx_dvec(len*len);
    for (int ii=0; ii < (len); ii++) {
        arma::vec temp_col = V.col(ii);
        for (int i=0; i < (len); i++) {
            std::complex<double> temp_im = im * ((dt*temp_col(i)) / 2.);
            a(pair_to_single(i, ii, len)) = 1. + 4.*r + temp_im;
            b(pair_to_single(i, ii, len)) = 1. - 4.*r - temp_im;
        }
    }
    return std::vector<arma::sp_cx_dmat>{create_mat(a, -r, len), create_mat(b, r, len)};
}

arma::sp_cx_dmat Diff::create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len) {
    int lenlen = len*len;
    arma::sp_cx_dmat A = arma::sp_cx_dmat(lenlen, lenlen);
    for (int i = 0; i < len; i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start, start, end, end) = create_tri(a, r, len, i);
    }
    for (int i = 1; i<(len); i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start-len, start, end-len, end) = create_rdiag(r, len);
        A.submat(start, start-len, end, end-len) = create_rdiag(r, len);
    }
    return A;
}

arma::sp_cx_dmat Diff::create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    for (int ii = 0; ii<len; ii++) {
        temp(ii, ii) = a(pair_to_single(ii, i, len));
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }
    return temp;
}

arma::sp_cx_dmat Diff::create_rdiag(const std::complex<double> r, const int len) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    temp.diag().fill(r);
    return temp;
}


Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    M = std::round(1./h);
    len = (M-2);
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    this->u = wave_init(xc, yc, widthx, widthy, px, py, M);
    AB = create_AB(V, h, dt, M);
    n_timesteps = std::round(T/dt);
    storage = arma::cube(len, len, n_timesteps + 1);
    Re = arma::cube(len, len, n_timesteps + 1);
    Im = arma::cube(len, len, n_timesteps + 1);
    this->T = T;
    this->dt = dt;
}

void Diff::run() {
    storage.slice(0) = this->probability(this->u);
    print_u(0);
    for (int i=1; i < n_timesteps + 1 ; i++) {
        b = AB.at(1) * this->u;
        this->u = arma::spsolve(AB.at(0), b);
        storage.slice(i) = this->probability(this->u);
        print_u(i);
    }
}

arma::mat Diff::probability(arma::cx_dvec &u) {
    arma::mat prob = arma::mat(len, len);
    for (int j=0; j < len; j++) {
        for (int i=0; i < len; i++) {
            prob(i, j) = std::real(std::conj(this->u(Diff::pair_to_single(i, j, len))) * this->u(Diff::pair_to_single(i, j, len)));
        }
    }
    return prob;
}

void Diff::print_u(int t) {

    arma::cx_mat temp = arma::cx_mat(len,len);
    int temp_len = len*len;
    for (int i = 0; i < temp_len; i++) {
        std::tuple<int,int> a = single_to_pair(i,len);
        temp(std::get<0>(a),std::get<1>(a)) = this->u(i);
        }
    Re.slice(t) = arma::real(temp);
    Im.slice(t) = arma::imag(temp);

}

arma::cx_dvec Diff::wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
    int len = M-2;
    arma::cx_dvec u = arma::cx_dvec(len*len);
    double h = 1./len;
    std::complex<double> im(0., 1.);
    for (int i=0; i<(len); i++) {
        for (int j=0; j<(len); j++) {
            double temp_x = ((h*(j+1)) - centerx);
            double temp_y = ((h*(i+1)) - centery);
            std::complex<double> arg = std::exp(-(temp_x*temp_x/(2*widthx*widthx)) - (temp_y*temp_y/(2*widthy*widthy)) + im*(px*temp_x) + im*(py*temp_y));
            u(pair_to_single(i, j, len)) = arg;
        }
    }
    u = u / std::sqrt(arma::accu(arma::conj(u) % u));

    return u;
}

arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M-2;
    double h = 1./len;
    arma::mat V = arma::mat(len, len).zeros();
    int j_start = std::round((centerx - (thickx/2))/h);
    int j_end = std::round((centerx + (thickx/2))/h) - 1;
    bool center_open = n_slit%2 != 0;   
    int counter = n_slit;
    int i_center = std::round(0.5/h) -1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture/h);
    int len_slit_sep = std::round(slit_sep/h);
    if (center_open) {
        int i_loc_up = i_center - std::round(aperture/(2*h));
        int i_loc_down = i_center + std::round(aperture/(2*h));
        i_loc_up--;
        i_loc_down++;
        counter--; 
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }
    else if (n_slit == 0) {
        for (int j=j_start; j<j_end; j++) {
            V.col(j).fill(potential);
        }
        return V;
    }
    else {
        int temp = std::round(len_slit_sep/2.);
        for (int ii=0; ii<temp; ii++) {
            for (int j=j_start; j<j_end; j++) {
                V(i_loc_up, j) = potential;
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }
}*/









/*Crank::Crank(double h, double dt, double T, double xc, double sigmax, double px,
             double yc, double sigmay, double py, double v0)
    : h(h), dt(dt), T(T), xc(xc), sigmax(sigmax), px(px), yc(yc), sigmay(sigmay),
      py(py), v0(v0) {
    U = arma::cx_cube(M, M, Nt + 1, arma::fill::zeros);
    V = arma::mat(M, M, arma::fill::zeros);
    A = arma::mat(M * M, M * M, arma::fill::zeros);
    B = arma::mat(M * M, M * M, arma::fill::zeros);
    Re = arma::cube(M, M, Nt + 1, arma::fill::zeros);
    Im = arma::cube(M, M, Nt + 1, arma::fill::zeros);

    // Set up potential, initial state, and matrices
    setupPotential();
    setupInitialState();
    setupMatrices();
}

// Function to set up the potential matrix V
void Crank::setupPotential() {
    // Initialize V to zeros
    V.zeros();

    // Parameters for the double-slit setup
    double wallThicknessX = 0.02;
    double wallPositionX = 0.5;
    double wallPieceLength = 0.05;
    double slitAperture = 0.05;

    // Loop over grid points and set potential based on the wall and slits
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            double x = i * h;
            double y = j * h;

            // Check if the point is within the wall
            if (x >= wallPositionX - 0.5 * wallThicknessX &&
                x <= wallPositionX + 0.5 * wallThicknessX) {
                // Check if the point is outside the slit region
                if (y < 0.5 - 0.5 * slitAperture || y > 0.5 + 0.5 * slitAperture) {
                    V(i, j) = std::numeric_limits<double>::infinity();
                }
            }
        }
    }
}

void Crank::setupInitialState() {
    // Initialize U0 to zeros
    U0.zeros();

    // Parameters for the initial wave packet
    double xc = 0.5;
    double yc = 0.5;
    double sigmaX = 0.05;
    double sigmaY = 0.05;
    double px = 1.0;
    double py = 0.0;

    // Loop over grid points and set initial state based on the Gaussian wave packet
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < M; ++j) {
            double x = i * h;
            double y = j * h;

            // Calculate Gaussian wave function
            U0(i, j) = exp(-(pow(x - xc, 2) / (2 * pow(sigmaX, 2))) -
                           (pow(y - yc, 2) / (2 * pow(sigmaY, 2)))) *
                        exp(cx(0, 1) * (px * x + py * y));
        }
    }
}

void Crank::setupMatrices() {
    // Calculate constants
    double r = dt / (2 * pow(h, 2));

    // Setup matrices A and B
    A.zeros();
    B.zeros();

    for (int k = 0; k < M - 2; ++k) {
        double a_k = 1.0 + 4.0 * r + cx(0, 1) * dt * V(k + 1, 0);
        double b_k = 1.0 - 4.0 * r - cx(0, 1) * dt * V(k + 1, 0);

        // Construct A_k and B_k matrices
        Mat<cx_double> A_k = a_k * eye<cx_mat>(M - 2, M - 2) - r * D(-1.0);
        Mat<cx_double> B_k = b_k * eye<cx_mat>(M - 2, M - 2) + r * D(-1.0);

        // Insert A_k and B_k into matrices A and B
        A.submat(span(k * (M - 2), (k + 1) * (M - 2) - 1),
                 span(k * (M - 2), (k + 1) * (M - 2) - 1)) = A_k;

        B.submat(span(k * (M - 2), (k + 1) * (M - 2) - 1),
                 span(k * (M - 2), (k + 1) * (M - 2) - 1)) = B_k;
    }
}


void Crank::runSimulation() {
    for (int n = 0; n < Nt; ++n) {
        Vec<cx_double> b = B * U;
        U = solve(A, b);
    }
}*/


/*int Diff::pair_to_single(const int i, const int j, const int len) {
    return (j) * len + (i);
}

std::tuple<int, int> Diff::single_to_pair(const int k, const int len) {
    int i = k % len;
    int j = k / len;
    return std::tuple<int, int>{i, j};
}

std::vector<arma::sp_cx_dmat> Diff::create_AB(arma::mat &V, const double h, const double dt, const int M) {
    int len = M-2;
    std::complex<double> im(0., 1.);
    std::complex<double> r = im*dt/(2*h*h);
    arma::cx_dvec a = arma::cx_dvec(len*len);
    arma::cx_dvec b = arma::cx_dvec(len*len);
    for (int ii=0; ii < (len); ii++) {
        arma::vec temp_col = V.col(ii);
        for (int i=0; i < (len); i++) {
            std::complex<double> temp_im = im * ((dt*temp_col(i)) / 2.);
            a(pair_to_single(i, ii, len)) = 1. + 4.*r + temp_im;
            b(pair_to_single(i, ii, len)) = 1. - 4.*r - temp_im;
        }
    }
    return std::vector<arma::sp_cx_dmat>{create_mat(a, -r, len), create_mat(b, r, len)};
}

arma::sp_cx_dmat Diff::create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len) {
    int lenlen = len*len;
    arma::sp_cx_dmat A = arma::sp_cx_dmat(lenlen, lenlen);
    for (int i = 0; i < len; i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start, start, end, end) = create_tri(a, r, len, i);
    }
    for (int i = 1; i<(len); i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start-len, start, end-len, end) = create_rdiag(r, len);
        A.submat(start, start-len, end, end-len) = create_rdiag(r, len);
    }
    return A;
}

arma::sp_cx_dmat Diff::create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    for (int ii = 0; ii<len; ii++) {
        temp(ii, ii) = a(pair_to_single(ii, i, len));
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }
    return temp;
}

arma::sp_cx_dmat Diff::create_rdiag(const std::complex<double> r, const int len) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    temp.diag().fill(r);
    return temp;
}


Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    M = std::round(1./h);
    len = (M-2);
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    this->u = wave_init(xc, yc, widthx, widthy, px, py, M);
    AB = create_AB(V, h, dt, M);
    n_timesteps = std::round(T/dt);
    storage = arma::cube(len, len, n_timesteps + 1);
    Re = arma::cube(len, len, n_timesteps + 1);
    Im = arma::cube(len, len, n_timesteps + 1);
    this->T = T;
    this->dt = dt;
}

void Diff::run() {
    storage.slice(0) = this->probability(this->u);
    print_u(0);
    for (int i=1; i < n_timesteps + 1 ; i++) {
        b = AB.at(1) * this->u;
        this->u = arma::spsolve(AB.at(0), b);
        storage.slice(i) = this->probability(this->u);
        print_u(i);
    }
}

arma::mat Diff::probability(arma::cx_dvec &u) {
    arma::mat prob = arma::mat(len, len);
    for (int j=0; j < len; j++) {
        for (int i=0; i < len; i++) {
            prob(i, j) = std::real(std::conj(this->u(Diff::pair_to_single(i, j, len))) * this->u(Diff::pair_to_single(i, j, len)));
        }
    }
    return prob;
}

void Diff::print_u(int t) {

    arma::cx_mat temp = arma::cx_mat(len,len);
    int temp_len = len*len;
    for (int i = 0; i < temp_len; i++) {
        std::tuple<int,int> a = single_to_pair(i,len);
        temp(std::get<0>(a),std::get<1>(a)) = this->u(i);
        }
    Re.slice(t) = arma::real(temp);
    Im.slice(t) = arma::imag(temp);

}

arma::cx_dvec Diff::wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
    int len = M-2;
    arma::cx_dvec u = arma::cx_dvec(len*len);
    double h = 1./len;
    std::complex<double> im(0., 1.);
    for (int i=0; i<(len); i++) {
        for (int j=0; j<(len); j++) {
            double temp_x = ((h*(j+1)) - centerx);
            double temp_y = ((h*(i+1)) - centery);
            std::complex<double> arg = std::exp(-(temp_x*temp_x/(2*widthx*widthx)) - (temp_y*temp_y/(2*widthy*widthy)) + im*(px*temp_x) + im*(py*temp_y));
            u(pair_to_single(i, j, len)) = arg;
        }
    }
    u = u / std::sqrt(arma::accu(arma::conj(u) % u));

    return u;
}

arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M-2;
    double h = 1./len;
    arma::mat V = arma::mat(len, len).zeros();
    int j_start = std::round((centerx - (thickx/2))/h);
    int j_end = std::round((centerx + (thickx/2))/h) - 1;
    bool center_open = n_slit%2 != 0;   
    int counter = n_slit;
    int i_center = std::round(0.5/h) -1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture/h);
    int len_slit_sep = std::round(slit_sep/h);
    if (center_open) {
        int i_loc_up = i_center - std::round(aperture/(2*h));
        int i_loc_down = i_center + std::round(aperture/(2*h));
        i_loc_up--;
        i_loc_down++;
        counter--; 
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }
    else if (n_slit == 0) {
        for (int j=j_start; j<j_end; j++) {
            V.col(j).fill(potential);
        }
        return V;
    }
    else {
        int temp = std::round(len_slit_sep/2.);
        for (int ii=0; ii<temp; ii++) {
            for (int j=j_start; j<j_end; j++) {
                V(i_loc_up, j) = potential;
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }
}*/



/*Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double width_x, double width_y, double potential, int n_slit, double thick_x, double center_x, double slit_sep, double aperture)
    : dt(dt), T(T), len(n_slit), M(2 * n_slit + 2), n_timesteps(static_cast<int>(T / dt)) {
    
    // Initialize cubes
    storage = arma::cube(M, M, n_timesteps, arma::fill::zeros);
    Re = arma::cube(M, M, n_timesteps, arma::fill::zeros);
    Im = arma::cube(M, M, n_timesteps, arma::fill::zeros);

    // Initialize spatial grid
    arma::vec x = arma::linspace(-len, len, M) * h;
    arma::mat X, Y;
    X = arma::repmat(x, 1, M);
    Y = arma::repmat(x.t(), M, 1);


    // Initialize initial state
    arma::cx_mat u_0(M, M, arma::fill::zeros);
    initializeInitialState(u_0, xc, yc, width_x, width_y, px, py, h);
    storage.slice(0) = real(u_0);
    Re.slice(0) = real(u_0);
    Im.slice(0) = imag(u_0);

    // Initialize potential
    arma::mat V(M, M, arma::fill::zeros);
    initializePotential(V, h, thick_x, center_x, slit_sep, aperture, potential);
}

inline int Diff::index2Dto1D(int i, int j, int M) {
    return i * M + j;
}

void Diff::generateMatrices(double r, arma::mat& A, arma::mat& B) {
    int vectorLength = (M - 2) * (M - 2);
    arma::vec a(vectorLength, arma::fill::zeros);
    arma::vec b(vectorLength, arma::fill::zeros);
    A.set_size(vectorLength, vectorLength);
    B.set_size(vectorLength, vectorLength);

    for (int i = 0; i < M - 2; ++i) {
        for (int j = 0; j < M - 2; ++j) {
            int k = index2Dto1D(i, j, M);
            a(k) = 1 + 4 * r;
            b(k) = 1 - 4 * r;
            if (i > 0) A(k, index2Dto1D(i - 1, j, M)) = -r;
            if (i < M - 3) A(k, index2Dto1D(i + 1, j, M)) = -r;
            A(k, k) = a(k);
            if (i > 0) B(k, index2Dto1D(i - 1, j, M)) = r;
            if (i < M - 3) B(k, index2Dto1D(i + 1, j, M)) = r;
            B(k, k) = b(k);
        }
    }
}

void Diff::crankNicolsonStep(arma::mat& A, arma::mat& B, arma::cx_mat& u_n) {
    int M = A.n_rows;
    arma::cx_mat b = B * u_n;
    for (int j = 0; j < M; ++j) {
        arma::cx_vec u_n_col = u_n.col(j);
        arma::cx_vec b_col = b.col(j);
        arma::cx_mat u_n1_col = arma::solve(A, arma::conv_to<arma::cx_mat>::from(b_col));



        u_n.col(j) = u_n1_col;
    }
}

void Diff::initializeInitialState(arma::cx_mat& u_0, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double h) {
    size_t M = u_0.n_rows;
    size_t N = u_0.n_cols;
    arma::vec x = arma::linspace<arma::vec>(0, (M - 1) * h, M);
    arma::vec y = arma::linspace<arma::vec>(0, (N - 1) * h, N);

    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            arma::cx_double exponent = -0.5 * ((x(i) - x_c) * (x(i) - x_c) / (sigma_x * sigma_x)
                + (y(j) - y_c) * (y(j) - y_c) / (sigma_y * sigma_y))
                + arma::cx_double(0, 1) * (p_x * x(i) + p_y * y(j));

            u_0(i, j) = std::exp(exponent);

            if (i == 0 || j == 0 || i == M - 1 || j == N - 1) {
                u_0(i, j) = arma::cx_double(0.0, 0.0);
            }
        }
    }

    double normalization = std::real(arma::accu(arma::conj(u_0) % u_0)) * h * h;
    u_0 /= std::sqrt(normalization);
}

void Diff::initializePotential(arma::mat& V, double h, double slitWidth, double wallThickness, double wallPosition, double slitSeparation, double barrierHeight) {
    int M = V.n_rows;
    int N = V.n_cols;
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            double x = i * h;
            double y = j * h;
            bool insideBarrier = (std::abs(x - wallPosition) < wallThickness / 2.0)
                && (std::abs(y) < slitSeparation / 2.0);
            V(i, j) = insideBarrier ? barrierHeight : 0.0;
        }
    }
}

void Diff::run() {
    for (int i = 1; i <= n_timesteps; ++i) {
        arma::mat A, B;
        generateMatrices(dt / (2 * h * h), A, B);
        crankNicolsonStep(A, B, storage.slice(i - 1));

        storage.slice(i) = real(Re.slice(i)) % real(Re.slice(i)) + imag(Im.slice(i)) % imag(Im.slice(i));
        Re.slice(i) = real(Re.slice(i));
        Im.slice(i) = imag(Im.slice(i));
    }
}*/


/*Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double width_x, double width_y, double potential, int n_slit, double thick_x, double center_x, double slit_sep, double aperture)
    : dt(dt), T(T), len(n_slit), M(2 * n_slit + 2), n_timesteps(static_cast<int>(T / dt)) {
    
    // Initialize cubes
    storage = arma::cube(M, M, n_timesteps, arma::fill::zeros);
    Re = arma::cube(M, M, n_timesteps, arma::fill::zeros);
    Im = arma::cube(M, M, n_timesteps, arma::fill::zeros);

    // Initialize spatial grid
    arma::vec x = arma::linspace(-len, len, M) * h;
    arma::mat X, Y;
    meshgrid(X, Y, x, x);

    // Initialize initial state
    arma::cx_mat u_0(M, M, arma::fill::zeros);
    initializeInitialState(u_0, xc, yc, width_x, width_y, px, py, h);
    storage.slice(0) = real(u_0);
    Re.slice(0) = real(u_0);
    Im.slice(0) = imag(u_0);

    // Initialize potential
    arma::mat V(M, M, arma::fill::zeros);
    initializePotential(V, h, thick_x, center_x, slit_sep, aperture, potential);
}

inline int Diff::index2Dto1D(int i, int j, int M) {
    return i * M + j;
}

void Diff::generateMatrices(int M, double r, arma::mat& A, arma::mat& B) {
    int vectorLength = (M - 2) * (M - 2);
    arma::vec a(vectorLength, arma::fill::zeros);
    arma::vec b(vectorLength, arma::fill::zeros);
    A.set_size(vectorLength, vectorLength);
    B.set_size(vectorLength, vectorLength);

    for (int i = 0; i < M - 2; ++i) {
        for (int j = 0; j < M - 2; ++j) {
            int k = index2Dto1D(i, j, M);
            a(k) = 1 + 4 * r;
            b(k) = 1 - 4 * r;
            if (i > 0) A(k, index2Dto1D(i - 1, j, M)) = -r;
            if (i < M - 3) A(k, index2Dto1D(i + 1, j, M)) = -r;
            A(k, k) = a(k);
            if (i > 0) B(k, index2Dto1D(i - 1, j, M)) = r;
            if (i < M - 3) B(k, index2Dto1D(i + 1, j, M)) = r;
            B(k, k) = b(k);
        }
    }
}

void Diff::crankNicolsonStep(arma::mat& A, arma::mat& B, arma::cx_mat& u_n) {
    int M = A.n_rows;
    arma::cx_mat b = B * u_n;
    for (int j = 0; j < M; ++j) {
        arma::cx_vec u_n_col = u_n.col(j);
        arma::cx_vec b_col = b.col(j);
        arma::cx_vec u_n1_col = arma::solve(A, b_col);
        u_n.col(j) = u_n1_col;
    }
}

void Diff::initializeInitialState(arma::cx_mat& u_0, double x_c, double y_c, double sigma_x, double sigma_y, double p_x, double p_y, double h) {
    size_t M = u_0.n_rows;
    size_t N = u_0.n_cols;
    arma::vec x = arma::linspace<arma::vec>(0, (M - 1) * h, M);
    arma::vec y = arma::linspace<arma::vec>(0, (N - 1) * h, N);

    for (size_t i = 0; i < M; ++i) {
        for (size_t j = 0; j < N; ++j) {
            arma::cx_double exponent = -0.5 * ((x(i) - x_c) * (x(i) - x_c) / (sigma_x * sigma_x)
                + (y(j) - y_c) * (y(j) - y_c) / (sigma_y * sigma_y))
                + arma::cx_double(0, 1) * (p_x * x(i) + p_y * y(j));

            u_0(i, j) = std::exp(exponent);

            if (i == 0 || j == 0 || i == M - 1 || j == N - 1) {
                u_0(i, j) = arma::cx_double(0.0, 0.0);
            }
        }
    }

    double normalization = std::real(arma::accu(arma::conj(u_0) % u_0)) * h * h;
    u_0 /= std::sqrt(normalization);
}

void Diff::initializePotential(arma::mat& V, double h, double slitWidth, double wallThickness, double wallPosition, double slitSeparation, double barrierHeight) {
    int M = V.n_rows;
    int N = V.n_cols;
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            double x = i * h;
            double y = j * h;
            bool insideBarrier = (std::abs(x - wallPosition) < wallThickness / 2.0)
                && (std::abs(y) < slitSeparation / 2.0);
            V(i, j) = insideBarrier ? barrierHeight : 0.0;
        }
    }
}

void Diff::run() {
    for (int i = 1; i <= n_timesteps; ++i) {
        arma::mat A, B;
        double r = dt / (2 * h * h);
        generateMatrices(M, r, A, B);

        crankNicolsonStep(A, B, storage.slice(i - 1));

        storage.slice(i) = Re.slice(i) % Re.slice(i) + Im.slice(i) % Im.slice(i);
        Re.slice(i) = real(storage.slice(i));
        Im.slice(i) = imag(storage.slice(i));
    }
}*/






/*
Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    M = std::round(1. / h);
    len = (M - 2);
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    this->u = wave_init(xc, yc, widthx, widthy, px, py, M);
    AB = create_AB(V, h, dt, M);
    n_timesteps = std::round(T / dt);
    storage = arma::cube(len, len, n_timesteps + 1);
    Re = arma::cube(len, len, n_timesteps + 1);
    Im = arma::cube(len, len, n_timesteps + 1);
    this->T = T;
    this->dt = dt;
}
int Diff::pair_to_single(const int i, const int j, const int len) {
    return (j) * len + (i);
}

std::tuple<int, int> Diff::single_to_pair(const int k, const int len) {
    int i = k % len;
    int j = k / len;
    return std::tuple<int, int>{i, j};
}
std::vector<arma::sp_cx_dmat> Diff::create_AB(arma::mat &V, const double h, const double dt, const int M) {
    int len = M-2;
    std::complex<double> im(0., 1.);
    std::complex<double> r = im*dt/(2*h*h);
    arma::cx_dvec a = arma::cx_dvec(len*len);
    arma::cx_dvec b = arma::cx_dvec(len*len);
    for (int ii=0; ii < (len); ii++) {
        arma::vec temp_col = V.col(ii);
        for (int i=0; i < (len); i++) {
            std::complex<double> temp_im = im * ((dt*temp_col(i)) / 2.);
            a(pair_to_single(i, ii, len)) = 1. + 4.*r + temp_im;
            b(pair_to_single(i, ii, len)) = 1. - 4.*r - temp_im;
        }
    }
    return std::vector<arma::sp_cx_dmat>{create_mat(a, -r, len), create_mat(b, r, len)};
}

arma::sp_cx_dmat Diff::create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len) {
    int lenlen = len * len;
    arma::sp_cx_dmat A = arma::sp_cx_dmat(lenlen, lenlen);

    for (int i = 0; i < len; i++) {
        int start = len * i;
        int end = len * i + (len - 1);
        A.submat(start, start, end, end) = create_tri(a, r, len, i);
    }

    for (int i = 1; i < len; i++) {
        int start = len * i;
        int end = len * i + (len - 1);
        A.submat(start - len, start, end - len, end) = create_rdiag(r, len);
        A.submat(start, start - len, end, end - len) = create_rdiag(r, len);
    }

    return A;
}
arma::sp_cx_dmat Diff::create_rdiag(const std::complex<double> r, const int len) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    temp.diag().fill(r);
    return temp;
}
void Diff::run() {
    storage.slice(0) = this->probability(this->u);
    print_u(0);
    for (int i = 1; i < n_timesteps + 1; i++) {
        b = AB.at(1) * this->u;
        this->u = arma::spsolve(AB.at(0), b);
        storage.slice(i) = this->probability(this->u);
        print_u(i);
    }
}
arma::sp_cx_dmat Diff::create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);

    for (int ii = 0; ii < len; ii++) {
        temp(ii, ii) = a(ii * len + i);
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }

    return temp;
}
arma::mat Diff::probability(arma::cx_dvec &u) {
    arma::mat prob = arma::mat(len, len);
    for (int j = 0; j < len; j++) {
        for (int i = 0; i < len; i++) {
            prob(i, j) = std::real(std::conj(this->u(j * len + i)) * this->u(j * len + i));
        }
    }
    return prob;
}

void Diff::print_u(int t) {
    arma::cx_mat temp = arma::cx_mat(len, len);
    int temp_len = len * len;
    for (int i = 0; i < temp_len; i++) {
        temp(i % len, i / len) = this->u(i);
    }
    Re.slice(t) = arma::real(temp);
    Im.slice(t) = arma::imag(temp);
}

arma::cx_dvec Diff::wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
    int len = M - 2;
    arma::cx_dvec u = arma::cx_dvec(len * len);
    double h = 1. / len;
    std::complex<double> im(0., 1.);
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            double temp_x = ((h * (j + 1)) - centerx);
            double temp_y = ((h * (i + 1)) - centery);
            std::complex<double> arg = std::exp(-(temp_x * temp_x / (2 * widthx * widthx)) -
                                                (temp_y * temp_y / (2 * widthy * widthy)) + im * (px * temp_x) + im * (py * temp_y));
            u(j * len + i) = arg;
        }
    }
    u = u / std::sqrt(arma::accu(arma::conj(u) % u));

    return u;
}

arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M - 2;
    double h = 1. / len;
    arma::mat V = arma::mat(len, len).zeros();
    int j_start = std::round((centerx - (thickx / 2)) / h);
    int j_end = std::round((centerx + (thickx / 2)) / h) - 1;
    bool center_open = n_slit % 2 != 0;

    int counter = n_slit;
    int i_center = std::round(0.5 / h) - 1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture / h);
    int len_slit_sep = std::round(slit_sep / h);

    if (center_open) {
        int i_loc_up = i_center - std::round(aperture / (2 * h));
        int i_loc_down = i_center + std::round(aperture / (2 * h));
        i_loc_up--;
        i_loc_down++;
        counter--;
        while (counter > 0) {
            for (int ii = 0; ii < len_slit_sep; ii++) {
                for (int j = j_start; j < j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i = i_loc_up; i >= 0; i--) {
            for (int j = j_start; j < j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i = i_loc_down; i < len; i++) {
            for (int j = j_start; j < j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    } else if (n_slit == 0) {
        for (int j = j_start; j < j_end; j++) {
            V.col(j).fill(potential);
        }
        return V;
    } else {
        int temp = std::round(len_slit_sep / 2.);
        for (int ii = 0; ii < temp; ii++) {
            for (int j = j_start; j < j_end; j++) {
                V(i_loc_up, j) = potential;
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;
        while (counter > 0) {
            for (int ii = 0; ii < len_slit_sep; ii++) {
                for (int j = j_start; j < j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i = i_loc_up; i >= 0; i--) {
            for (int j = j_start; j < j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i = i_loc_down; i < len; i++) {
            for (int j = j_start; j < j_end; j++) {
                V(i, j) = potential;
            }
        }

        return V;
    }
}*/










///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def
///////////////////////////Two def

/*int Diff::pair_to_single(const int i, const int j, const int len) {
    return (j) * len + (i);
}

std::tuple<int, int> Diff::single_to_pair(const int k, const int len) {
    int i = k % len;
    int j = k / len;
    return std::tuple<int, int>{i, j};
}

std::vector<arma::sp_cx_dmat> Diff::create_AB(arma::mat &V, const double h, const double dt, const int M) {
    int len = M-2;
    std::complex<double> im(0., 1.);
    std::complex<double> r = im*dt/(2*h*h);
    arma::cx_dvec a = arma::cx_dvec(len*len);
    arma::cx_dvec b = arma::cx_dvec(len*len);
    for (int ii=0; ii < (len); ii++) {
        arma::vec temp_col = V.col(ii);
        for (int i=0; i < (len); i++) {
            std::complex<double> temp_im = im * ((dt*temp_col(i)) / 2.);
            a(pair_to_single(i, ii, len)) = 1. + 4.*r + temp_im;
            b(pair_to_single(i, ii, len)) = 1. - 4.*r - temp_im;
        }
    }
    return std::vector<arma::sp_cx_dmat>{create_mat(a, -r, len), create_mat(b, r, len)};
}

arma::sp_cx_dmat Diff::create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len) {
    int lenlen = len*len;
    arma::sp_cx_dmat A = arma::sp_cx_dmat(lenlen, lenlen);
    for (int i = 0; i < len; i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start, start, end, end) = create_tri(a, r, len, i);
    }
    for (int i = 1; i<(len); i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start-len, start, end-len, end) = create_rdiag(r, len);
        A.submat(start, start-len, end, end-len) = create_rdiag(r, len);
    }
    return A;
}
arma::sp_cx_dmat Diff::create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    for (int ii = 0; ii<len; ii++) {
        temp(ii, ii) = a(pair_to_single(ii, i, len));
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }
    return temp;
}

arma::sp_cx_dmat Diff::create_rdiag(const std::complex<double> r, const int len) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    temp.diag().fill(r);
    return temp;
}


Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    M = std::round(1./h);
    len = (M-2);
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    this->u = wave_init(xc, yc, widthx, widthy, px, py, M);
    AB = create_AB(V, h, dt, M);
    n_timesteps = std::round(T/dt);
    storage = arma::cube(len, len, n_timesteps + 1);
    Re = arma::cube(len, len, n_timesteps + 1);
    Im = arma::cube(len, len, n_timesteps + 1);
    this->T = T;
    this->dt = dt;
}

void Diff::run() {
    storage.slice(0) = this->probability(this->u);
    print_u(0);
    for (int i=1; i < n_timesteps + 1 ; i++) {
        b = AB.at(1) * this->u;
        this->u = arma::spsolve(AB.at(0), b);
        storage.slice(i) = this->probability(this->u);
        print_u(i);
    }
}

arma::mat Diff::probability(arma::cx_dvec &u) {
    arma::mat prob = arma::mat(len, len);
    for (int j=0; j < len; j++) {
        for (int i=0; i < len; i++) {
            prob(i, j) = std::real(std::conj(this->u(Diff::pair_to_single(i, j, len))) * this->u(Diff::pair_to_single(i, j, len)));
        }
    }
    return prob;
}

void Diff::print_u(int t) {

    arma::cx_mat temp = arma::cx_mat(len,len);
    int temp_len = len*len;
    for (int i = 0; i < temp_len; i++) {
        std::tuple<int,int> a = single_to_pair(i,len);
        temp(std::get<0>(a),std::get<1>(a)) = this->u(i);
        }
    Re.slice(t) = arma::real(temp);
    Im.slice(t) = arma::imag(temp);

}


arma::cx_dvec Diff::wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
    int len = M-2;
    arma::cx_dvec u = arma::cx_dvec(len*len);
    double h = 1./len;
    std::complex<double> im(0., 1.);
    for (int i=0; i<(len); i++) {
        for (int j=0; j<(len); j++) {
            double temp_x = ((h*(j+1)) - centerx);
            double temp_y = ((h*(i+1)) - centery);
            std::complex<double> arg = std::exp(-(temp_x*temp_x/(2*widthx*widthx)) - (temp_y*temp_y/(2*widthy*widthy)) + im*(px*temp_x) + im*(py*temp_y));
            u(pair_to_single(i, j, len)) = arg;
        }
    }
    u = u / std::sqrt(arma::accu(arma::conj(u) % u));

    return u;
}


arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M-2;
    double h = 1./len;
    arma::mat V = arma::mat(len, len).zeros();
    int j_start = std::round((centerx - (thickx/2))/h);
    int j_end = std::round((centerx + (thickx/2))/h) - 1;
    bool center_open = n_slit%2 != 0;   

    int counter = n_slit;
    int i_center = std::round(0.5/h) -1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture/h);
    int len_slit_sep = std::round(slit_sep/h);

    if (center_open) {
        int i_loc_up = i_center - std::round(aperture/(2*h));
        int i_loc_down = i_center + std::round(aperture/(2*h));
        i_loc_up--;
        i_loc_down++;
        counter--; 
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }

    else if (n_slit == 0) {
        for (int j=j_start; j<j_end; j++) {
            V.col(j).fill(potential);
        }
        return V;
    }
    else {
        int temp = std::round(len_slit_sep/2.);
        for (int ii=0; ii<temp; ii++) {
            for (int j=j_start; j<j_end; j++) {
                V(i_loc_up, j) = potential;
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }

        return V;
    }
}*/
































///////////////////////////Matrix def
///////////////////////////Matrix def
///////////////////////////Matrix def
///////////////////////////Matrix def
///////////////////////////Matrix def
///////////////////////////Matrix def///////////////////////////Matrix def
///////////////////////////Matrix def
///////////////////////////Matrix def
/*Diff::Diff(double h_in, double dt_in, double T_in)
{
  h = h_in;
  dt = dt_in;
  T = T_in;

  M = 1. / h + 1;
  N = T / dt + 1;
  r = arma::cx_double(0, dt / (2 * h * h));

  U.zeros(M - 2, M - 2);
  V.zeros(M - 2, M - 2);
  S.zeros(M - 2, M - 2, N);

  A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
  B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
}

// Get index k from (i, j)
int Diff::index_k(int i, int j)
{
  return i + (M - 2) * j;
}

// Set potential barriers
void Diff::set_potential(std::string filename)
{
  // Load the file into the potential
  V.load(filename, arma::raw_ascii);
  assert(V.n_rows == (M - 2) && V.n_cols == (M - 2));

}

// Fill the matrices according to the Crank-Nicolson regime
void Diff::fill_matrices()
{
  // Vectors a and b, which will be the main diagonals of A and B
  arma::cx_vec a((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_vec b((M - 2) * (M - 2), arma::fill::zeros);

  // Tridiagonal sub matrix with signature (r, 0, r)
  arma::cx_mat sub_diag(M - 2, M - 2, arma::fill::zeros);
  sub_diag.diag(-1).fill(r);
  sub_diag.diag(1).fill(r);

  // Filling diagonals of A and B
  A.diag(M - 2).fill(-r);
  A.diag(2 - M).fill(-r);
  B.diag(M - 2).fill(r);
  B.diag(2 - M).fill(r);

  // Filling a and b
  arma::cx_double a_k, b_k;
  int k;

  for (int i = 0; i < M - 2; i++)
  {
    for (int j = 0; j < M - 2; j++)
    {
      k = index_k(i, j);

      a_k = arma::cx_double(1., 4. * r.imag() + .5 * dt * V(i, j));
      b_k = std::conj(a_k);

      a(k) = a_k;
      b(k) = b_k;
    }
  }

  // Filling A and B with sub matrices
  for (int i = 0; i < M - 2; i++)
  {
    int j = i * (M - 2);

    A.submat(j, j, j + M - 3, j + M - 3) = -sub_diag;
    B.submat(j, j, j + M - 3, j + M - 3) = sub_diag;
  }

  // Setting main diagonals of A and B
  A.diag() = a;
  B.diag() = b;
}

// Solve matrix eq. Au^(n+1) = Bu^n
void Diff::solve()
{
  // Solving the equation for each time step and saving as slice in S
  for (int i = 1; i < N; i++)
  {
    u_current = U.as_col();
    u_new = arma::spsolve(A, B * u_current);

    // Refilling u vector into U matrix
    for (int n = 0; n < U.n_cols; n ++)
    {
      for (int k = 0; k < M - 2; k++)
      {
        U.col(n)(k) = u_new(k + n * (M - 2));
      }
    }

    S.slice(i) = U;
    u_current = u_new;
  }
}

// Set the initial state of the system
void Diff::set_initial_state(double x_c, double sigma_x, double p_x, double y_c,
                       double sigma_y, double p_y)
{
  arma::cx_double exponent;
  double re, im;
  int xi, yi;

  for (int i = 0; i < M - 2; i++)
  {
    yi = i;

    for (int j = 0; j < M - 2; j++)
    {
      xi = j;

      re = -(xi * h - x_c) * (xi * h * x_c) / (2 * sigma_x * sigma_x)
           -(yi * h - y_c) * (yi * h - y_c) / (2 * sigma_y * sigma_y);

      im = p_x * (xi * h - x_c) + (yi * h - y_c);

      exponent = arma::cx_double(re, im);

      U(i, j) = std::exp(exponent);
    }
  }

  // Normalize
  U /= std::sqrt(arma::accu(arma::conj(U) % U));
  S.slice(0) = U;
}*/



/*Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    M = 1. / h + 1;
    r = arma::cx_double(0, dt / (2 * h * h));
    //len = (M-2);
    //V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    //this->u = wave_init(xc, yc, widthx, widthy, px, py, M);
    //AB = create_AB(V, h, dt, M);
    //n_timesteps = std::round(T/dt);
    //storage = arma::cube(len, len, n_timesteps + 1);
    //Re = arma::cube(len, len, n_timesteps + 1);
    //Im = arma::cube(len, len, n_timesteps + 1);
    this->T = T;
    this->dt = dt;
    A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
}*/

/*arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M-2;
    double h = 1./len;
    arma::mat V = arma::mat(len, len).zeros();
    int j_start = std::round((centerx - (thickx/2))/h);
    int j_end = std::round((centerx + (thickx/2))/h) - 1;
    bool center_open = n_slit%2 != 0;   
    int counter = n_slit;
    int i_center = std::round(0.5/h) -1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture/h);
    int len_slit_sep = std::round(slit_sep/h);
    if (center_open) {
        int i_loc_up = i_center - std::round(aperture/(2*h));
        int i_loc_down = i_center + std::round(aperture/(2*h));
        i_loc_up--;
        i_loc_down++;
        counter--; 
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }
    else if (n_slit == 0) {
        for (int j=j_start; j<j_end; j++) {
            V.col(j).fill(potential);
        }
        return V;
    }
    else {
        int temp = std::round(len_slit_sep/2.);
        for (int ii=0; ii<temp; ii++) {
            for (int j=j_start; j<j_end; j++) {
                V(i_loc_up, j) = potential;
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }
}*/

/*int Diff::index_k(int i, int j){
  return i + (M - 2) * j;
}

void Diff::fill_matrices()
{
  // Vectors a and b, which will be the main diagonals of A and B
  arma::cx_vec a((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_vec b((M - 2) * (M - 2), arma::fill::zeros);

  // Tridiagonal sub matrix with signature (r, 0, r)
  arma::cx_mat sub_diag(M - 2, M - 2, arma::fill::zeros);
  sub_diag.diag(-1).fill(r);
  sub_diag.diag(1).fill(r);

  // Filling diagonals of A and B
  A.diag(M - 2).fill(-r);
  A.diag(2 - M).fill(-r);
  B.diag(M - 2).fill(r);
  B.diag(2 - M).fill(r);

  // Filling a and b
  arma::cx_double a_k, b_k;
  int k;

  for (int i = 0; i < M - 2; i++)
  {
    for (int j = 0; j < M - 2; j++)
    {
      k = index_k(i, j);

      a_k = arma::cx_double(1., 4. * r.imag() + .5 * dt * V(i, j));
      b_k = std::conj(a_k);

      a(k) = a_k;
      b(k) = b_k;
    }
  }

  // Filling A and B with sub matrices
  for (int i = 0; i < M - 2; i++)
  {
    int j = i * (M - 2);

    std::cout << "Submatrices dimensions: " << A.submat(j, j, j + M - 3, j + M - 3).n_rows << " x " << A.submat(j, j, j + M - 3, j + M - 3).n_cols << std::endl;

    A.submat(j, j, j + M - 3, j + M - 3) = -sub_diag;
    B.submat(j, j, j + M - 3, j + M - 3) = sub_diag;
  }

  // Setting main diagonals of A and B
  A.diag() = a;
  B.diag() = b;
}
*/



/*int Diff::index_k(int i, int j){
  return i + (M - 2) * j;
}
void Diff::fill_matrices(){
  arma::cx_vec a((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_vec b((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_mat sub_diag(M - 2, M - 2, arma::fill::zeros);
  sub_diag.diag(-1).fill(r);
  sub_diag.diag(1).fill(r);
  A.diag(M - 2).fill(-r);
  A.diag(2 - M).fill(-r);
  B.diag(M - 2).fill(r);
  B.diag(2 - M).fill(r);
  arma::cx_double a_k, b_k;
  int k;
  for (int i = 0; i < M - 2; i++){
    for (int j = 0; j < M - 2; j++){
      k = index_k(i, j);
      a_k = arma::cx_double(1., 4. * r.imag() + .5 * dt * V(i, j));
      b_k = std::conj(a_k);
      a(k) = a_k;
      b(k) = b_k;
    }
  }
  for (int i = 0; i < M - 2; i++){
    int j = i * (M - 2);
    A.submat(j, j, j + M - 3, j + M - 3) = -sub_diag;
    B.submat(j, j, j + M - 3, j + M - 3) = sub_diag;
  }
  A.diag() = a;
  B.diag() = b;
}*/





/*int Diff::pair_to_single(const int i, const int j, const int len) {
    return (j) * len + (i);
}

std::tuple<int, int> Diff::single_to_pair(const int k, const int len) {
    int i = k % len;
    int j = k / len;
    return std::tuple<int, int>{i, j};
}

std::vector<arma::sp_cx_dmat> Diff::create_AB(arma::mat &V, const double h, const double dt, const int M) {
    int len = M-2;
    std::complex<double> im(0., 1.);
    std::complex<double> r = im*dt/(2*h*h);
    arma::cx_dvec a = arma::cx_dvec(len*len);
    arma::cx_dvec b = arma::cx_dvec(len*len);
    for (int ii=0; ii < (len); ii++) {
        arma::vec temp_col = V.col(ii);
        for (int i=0; i < (len); i++) {
            std::complex<double> temp_im = im * ((dt*temp_col(i)) / 2.);
            a(pair_to_single(i, ii, len)) = 1. + 4.*r + temp_im;
            b(pair_to_single(i, ii, len)) = 1. - 4.*r - temp_im;
        }
    }
    return std::vector<arma::sp_cx_dmat>{create_mat(a, -r, len), create_mat(b, r, len)};
}

arma::sp_cx_dmat Diff::create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len) {
    int lenlen = len*len;
    arma::sp_cx_dmat A = arma::sp_cx_dmat(lenlen, lenlen);
    for (int i = 0; i < len; i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start, start, end, end) = create_tri(a, r, len, i);
    }
    for (int i = 1; i<(len); i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start-len, start, end-len, end) = create_rdiag(r, len);
        A.submat(start, start-len, end, end-len) = create_rdiag(r, len);
    }
    return A;
}
arma::sp_cx_dmat Diff::create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    for (int ii = 0; ii<len; ii++) {
        temp(ii, ii) = a(pair_to_single(ii, i, len));
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }
    return temp;
}

arma::sp_cx_dmat Diff::create_rdiag(const std::complex<double> r, const int len) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    temp.diag().fill(r);
    return temp;
}


Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    M = std::round(1./h);
    len = (M-2);
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    this->u = wave_init(xc, yc, widthx, widthy, px, py, M);
    AB = create_AB(V, h, dt, M);
    n_timesteps = std::round(T/dt);
    storage = arma::cube(len, len, n_timesteps + 1);
    Re = arma::cube(len, len, n_timesteps + 1);
    Im = arma::cube(len, len, n_timesteps + 1);
    this->T = T;
    this->dt = dt;

    M = 1. / h + 1;
    N = T / dt + 1;
    r = arma::cx_double(0, dt / (2 * h * h));

    U.zeros(M - 2, M - 2);
    V.zeros(M - 2, M - 2);
    S.zeros(M - 2, M - 2, N);

    A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    }

void Diff::run() {
    storage.slice(0) = this->probability(this->u);
    print_u(0);
    for (int i=1; i < n_timesteps + 1 ; i++) {
        b = AB.at(1) * this->u;
        this->u = arma::spsolve(AB.at(0), b);
        storage.slice(i) = this->probability(this->u);
        print_u(i);
    }
}



arma::mat Diff::probability(arma::cx_dvec &u) {
    arma::mat prob = arma::mat(len, len);
    for (int j=0; j < len; j++) {
        for (int i=0; i < len; i++) {
            prob(i, j) = std::real(std::conj(this->u(Diff::pair_to_single(i, j, len))) * this->u(Diff::pair_to_single(i, j, len)));
        }
    }
    return prob;
}

void Diff::print_u(int t) {

    arma::cx_mat temp = arma::cx_mat(len,len);
    int temp_len = len*len;
    for (int i = 0; i < temp_len; i++) {
        std::tuple<int,int> a = single_to_pair(i,len);
        temp(std::get<0>(a),std::get<1>(a)) = this->u(i);
        }
    Re.slice(t) = arma::real(temp);
    Im.slice(t) = arma::imag(temp);

}


arma::cx_dvec Diff::wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
    int len = M-2;
    arma::cx_dvec u = arma::cx_dvec(len*len);
    double h = 1./len;
    std::complex<double> im(0., 1.);
    for (int i=0; i<(len); i++) {
        for (int j=0; j<(len); j++) {
            double temp_x = ((h*(j+1)) - centerx);
            double temp_y = ((h*(i+1)) - centery);
            std::complex<double> arg = std::exp(-(temp_x*temp_x/(2*widthx*widthx)) - (temp_y*temp_y/(2*widthy*widthy)) + im*(px*temp_x) + im*(py*temp_y));
            u(pair_to_single(i, j, len)) = arg;
        }
    }
    u = u / std::sqrt(arma::accu(arma::conj(u) % u));

    return u;
}


arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M-2;
    double h = 1./len;
    arma::mat V = arma::mat(len, len).zeros();
    int j_start = std::round((centerx - (thickx/2))/h);
    int j_end = std::round((centerx + (thickx/2))/h) - 1;
    bool center_open = n_slit%2 != 0;   

    int counter = n_slit;
    int i_center = std::round(0.5/h) -1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture/h);
    int len_slit_sep = std::round(slit_sep/h);

    if (center_open) {
        int i_loc_up = i_center - std::round(aperture/(2*h));
        int i_loc_down = i_center + std::round(aperture/(2*h));
        i_loc_up--;
        i_loc_down++;
        counter--; 
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }

    else if (n_slit == 0) {
        for (int j=j_start; j<j_end; j++) {
            V.col(j).fill(potential);
        }
        return V;
    }
    else {
        int temp = std::round(len_slit_sep/2.);
        for (int ii=0; ii<temp; ii++) {
            for (int j=j_start; j<j_end; j++) {
                V(i_loc_up, j) = potential;
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }

        return V;
    }
}


int Diff::index_k(int i, int j){
  return i + (M - 2) * j;
}

void Diff::fill_matrices(){
  arma::cx_vec a((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_vec b((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_mat sub_diag(M - 2, M - 2, arma::fill::zeros);
  sub_diag.diag(-1).fill(r);
  sub_diag.diag(1).fill(r);
  A.diag(M - 2).fill(-r);
  A.diag(2 - M).fill(-r);
  B.diag(M - 2).fill(r);
  B.diag(2 - M).fill(r);
  arma::cx_double a_k, b_k;
  int k;
  for (int i = 0; i < M - 2; i++){
    for (int j = 0; j < M - 2; j++){
      k = index_k(i, j);
      a_k = arma::cx_double(1., 4. * r.imag() + .5 * dt * V(i, j));
      b_k = std::conj(a_k);
      a(k) = a_k;
      b(k) = b_k;
    }
  }
  for (int i = 0; i < M - 2; i++){
    int j = i * (M - 2);
    A.submat(j, j, j + M - 3, j + M - 3) = -sub_diag;
    B.submat(j, j, j + M - 3, j + M - 3) = sub_diag;
  }
  A.diag() = a;
  B.diag() = b;
}

void Diff::print_u(int t) {

    arma::cx_mat temp = arma::cx_mat(len,len);
    int temp_len = len*len;
    for (int i = 0; i < temp_len; i++) {
        std::tuple<int,int> a = single_to_pair(i,len);
        temp(std::get<0>(a),std::get<1>(a)) = this->u(i);
        }
    Re.slice(t) = arma::real(temp);
    Im.slice(t) = arma::imag(temp);

}

void Diff::solve(){
  for (int i = 1; i < N; i++){
    u_current = U.as_col();
    u_new = arma::spsolve(A, B * u_current);
    for (int n = 0; n < U.n_cols; n ++){
      for (int k = 0; k < M - 2; k++){
        U.col(n)(k) = u_new(k + n * (M - 2));
      }
    }
    print_u(i);
    S.slice(i) = U;
    u_current = u_new;
  }
}

void Diff::set_initial_state(double x_c, double sigma_x, double p_x, double y_c,double sigma_y, double p_y){
  arma::cx_double exponent;
  double re, im;
  int xi, yi;
  for (int i = 0; i < M - 2; i++){
    yi = i;
    for (int j = 0; j < M - 2; j++){
      xi = j;
      re = -(xi * h - x_c) * (xi * h * x_c) / (2 * sigma_x * sigma_x)-(yi * h - y_c) * (yi * h - y_c) / (2 * sigma_y * sigma_y);
      im = p_x * (xi * h - x_c) + (yi * h - y_c);
      exponent = arma::cx_double(re, im);
      U(i, j) = std::exp(exponent);

    }
  }
  U /= std::sqrt(arma::accu(arma::conj(U) % U));
  S.slice(0) = U;
}


int Diff::pair_to_single(const int i, const int j, const int len) {
    return (j) * len + (i);
}

std::tuple<int, int> Diff::single_to_pair(const int k, const int len) {
    int i = k % len;
    int j = k / len;
    return std::tuple<int, int>{i, j};
}

std::vector<arma::sp_cx_dmat> Diff::create_AB(arma::mat &V, const double h, const double dt, const int M) {
    int len = M-2;
    std::complex<double> im(0., 1.);
    std::complex<double> r = im*dt/(2*h*h);
    arma::cx_dvec a = arma::cx_dvec(len*len);
    arma::cx_dvec b = arma::cx_dvec(len*len);
    for (int ii=0; ii < (len); ii++) {
        arma::vec temp_col = V.col(ii);
        for (int i=0; i < (len); i++) {
            std::complex<double> temp_im = im * ((dt*temp_col(i)) / 2.);
            a(pair_to_single(i, ii, len)) = 1. + 4.*r + temp_im;
            b(pair_to_single(i, ii, len)) = 1. - 4.*r - temp_im;
        }
    }
    return std::vector<arma::sp_cx_dmat>{create_mat(a, -r, len), create_mat(b, r, len)};
}

arma::sp_cx_dmat Diff::create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len) {
    int lenlen = len*len;
    arma::sp_cx_dmat A = arma::sp_cx_dmat(lenlen, lenlen);
    for (int i = 0; i < len; i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start, start, end, end) = create_tri(a, r, len, i);
    }
    for (int i = 1; i<(len); i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start-len, start, end-len, end) = create_rdiag(r, len);
        A.submat(start, start-len, end, end-len) = create_rdiag(r, len);
    }
    return A;
}
arma::sp_cx_dmat Diff::create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    for (int ii = 0; ii<len; ii++) {
        temp(ii, ii) = a(pair_to_single(ii, i, len));
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }
    return temp;
}

arma::sp_cx_dmat Diff::create_rdiag(const std::complex<double> r, const int len) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    temp.diag().fill(r);
    return temp;
}


Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    M = std::round(1./h);
    len = (M-2);
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    this->u = wave_init(xc, yc, widthx, widthy, px, py, M);
    AB = create_AB(V, h, dt, M);
    n_timesteps = std::round(T/dt);
    storage = arma::cube(len, len, n_timesteps + 1);
    Re = arma::cube(len, len, n_timesteps + 1);
    Im = arma::cube(len, len, n_timesteps + 1);
    this->T = T;
    this->dt = dt;



    M = 1. / h + 1;
    N = T / dt + 1;
    r = arma::cx_double(0, dt / (2 * h * h));

    U.zeros(M - 2, M - 2);
    V.zeros(M - 2, M - 2);
    S.zeros(M - 2, M - 2, N);

    A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
    }

void Diff::run() {
    storage.slice(0) = this->probability(this->u);
    print_u(0);
    for (int i=1; i < n_timesteps + 1 ; i++) {
        b = AB.at(1) * this->u;
        this->u = arma::spsolve(AB.at(0), b);
        storage.slice(i) = this->probability(this->u);
        print_u(i);
    }
}



arma::mat Diff::probability(arma::cx_dvec &u) {
    arma::mat prob = arma::mat(len, len);
    for (int j=0; j < len; j++) {
        for (int i=0; i < len; i++) {
            prob(i, j) = std::real(std::conj(this->u(Diff::pair_to_single(i, j, len))) * this->u(Diff::pair_to_single(i, j, len)));
        }
    }
    return prob;
}

void Diff::print_u(int t) {

    arma::cx_mat temp = arma::cx_mat(len,len);
    int temp_len = len*len;
    for (int i = 0; i < temp_len; i++) {
        std::tuple<int,int> a = single_to_pair(i,len);
        temp(std::get<0>(a),std::get<1>(a)) = this->u(i);
        }
    Re.slice(t) = arma::real(temp);
    Im.slice(t) = arma::imag(temp);

}


arma::cx_dvec Diff::wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
    int len = M-2;
    arma::cx_dvec u = arma::cx_dvec(len*len);
    double h = 1./len;
    std::complex<double> im(0., 1.);
    for (int i=0; i<(len); i++) {
        for (int j=0; j<(len); j++) {
            double temp_x = ((h*(j+1)) - centerx);
            double temp_y = ((h*(i+1)) - centery);
            std::complex<double> arg = std::exp(-(temp_x*temp_x/(2*widthx*widthx)) - (temp_y*temp_y/(2*widthy*widthy)) + im*(px*temp_x) + im*(py*temp_y));
            u(pair_to_single(i, j, len)) = arg;
        }
    }
    u = u / std::sqrt(arma::accu(arma::conj(u) % u));

    return u;
}


arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M-2;
    double h = 1./len;
    arma::mat V = arma::mat(len, len).zeros();
    int j_start = std::round((centerx - (thickx/2))/h);
    int j_end = std::round((centerx + (thickx/2))/h) - 1;
    bool center_open = n_slit%2 != 0;   

    int counter = n_slit;
    int i_center = std::round(0.5/h) -1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture/h);
    int len_slit_sep = std::round(slit_sep/h);

    if (center_open) {
        int i_loc_up = i_center - std::round(aperture/(2*h));
        int i_loc_down = i_center + std::round(aperture/(2*h));
        i_loc_up--;
        i_loc_down++;
        counter--; 
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }

    else if (n_slit == 0) {
        for (int j=j_start; j<j_end; j++) {
            V.col(j).fill(potential);
        }
        return V;
    }
    else {
        int temp = std::round(len_slit_sep/2.);
        for (int ii=0; ii<temp; ii++) {
            for (int j=j_start; j<j_end; j++) {
                V(i_loc_up, j) = potential;
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }

        return V;
    }
}



int Diff::index_k(int i, int j)
{
  return i + (M - 2) * j;
}


void Diff::fill_matrices()
{
  // Vectors a and b, which will be the main diagonals of A and B
  arma::cx_vec a((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_vec b((M - 2) * (M - 2), arma::fill::zeros);

  // Tridiagonal sub matrix with signature (r, 0, r)
  arma::cx_mat sub_diag(M - 2, M - 2, arma::fill::zeros);
  sub_diag.diag(-1).fill(r);
  sub_diag.diag(1).fill(r);

  // Filling diagonals of A and B
  A.diag(M - 2).fill(-r);
  A.diag(2 - M).fill(-r);
  B.diag(M - 2).fill(r);
  B.diag(2 - M).fill(r);

  // Filling a and b
  arma::cx_double a_k, b_k;
  int k;

  for (int i = 0; i < M - 2; i++)
  {
    for (int j = 0; j < M - 2; j++)
    {
      k = index_k(i, j);

      a_k = arma::cx_double(1., 4. * r.imag() + .5 * dt * V(i, j));
      b_k = std::conj(a_k);

      a(k) = a_k;
      b(k) = b_k;
    }
  }

  // Filling A and B with sub matrices
  for (int i = 0; i < M - 2; i++)
  {
    int j = i * (M - 2);

    A.submat(j, j, j + M - 3, j + M - 3) = -sub_diag;
    B.submat(j, j, j + M - 3, j + M - 3) = sub_diag;
  }

  // Setting main diagonals of A and B
  A.diag() = a;
  B.diag() = b;
}

// Solve matrix eq. Au^(n+1) = Bu^n
void Diff::solve()
{
  // Solving the equation for each time step and saving as slice in S
  for (int i = 1; i < N; i++)
  {
    u_current = U.as_col();
    u_new = arma::spsolve(A, B * u_current);

    // Refilling u vector into U matrix
    for (int n = 0; n < U.n_cols; n ++)
    {
      for (int k = 0; k < M - 2; k++)
      {
        U.col(n)(k) = u_new(k + n * (M - 2));
      }
    }

    S.slice(i) = U;
    u_current = u_new;
  }
}

// Set the initial state of the system
void Diff::set_initial_state(double x_c, double sigma_x, double p_x, double y_c,
                       double sigma_y, double p_y)
{
  arma::cx_double exponent;
  double re, im;
  int xi, yi;

  for (int i = 0; i < M - 2; i++)
  {
    yi = i;

    for (int j = 0; j < M - 2; j++)
    {
      xi = j;

      re = -(xi * h - x_c) * (xi * h * x_c) / (2 * sigma_x * sigma_x)
           -(yi * h - y_c) * (yi * h - y_c) / (2 * sigma_y * sigma_y);

      im = p_x * (xi * h - x_c) + (yi * h - y_c);

      exponent = arma::cx_double(re, im);

      U(i, j) = std::exp(exponent);
    }
  }

  // Normalize
  U /= std::sqrt(arma::accu(arma::conj(U) % U));
  S.slice(0) = U;
}*/


/*int Diff::pair_to_single(const int i, const int j, const int len) {
    return (j) * len + (i);
}

std::tuple<int, int> Diff::single_to_pair(const int k, const int len) {
    int i = k % len;
    int j = k / len;
    return std::tuple<int, int>{i, j};
}

std::vector<arma::sp_cx_dmat> Diff::create_AB(arma::mat &V, const double h, const double dt, const int M) {
    int len = M-2;
    std::complex<double> im(0., 1.);
    std::complex<double> r = im*dt/(2*h*h);
    arma::cx_dvec a = arma::cx_dvec(len*len);
    arma::cx_dvec b = arma::cx_dvec(len*len);
    for (int ii=0; ii < (len); ii++) {
        arma::vec temp_col = V.col(ii);
        for (int i=0; i < (len); i++) {
            std::complex<double> temp_im = im * ((dt*temp_col(i)) / 2.);
            a(pair_to_single(i, ii, len)) = 1. + 4.*r + temp_im;
            b(pair_to_single(i, ii, len)) = 1. - 4.*r - temp_im;
        }
    }
    return std::vector<arma::sp_cx_dmat>{create_mat(a, -r, len), create_mat(b, r, len)};
}

arma::sp_cx_dmat Diff::create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len) {
    int lenlen = len*len;
    arma::sp_cx_dmat A = arma::sp_cx_dmat(lenlen, lenlen);
    for (int i = 0; i < len; i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start, start, end, end) = create_tri(a, r, len, i);
    }
    for (int i = 1; i<(len); i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start-len, start, end-len, end) = create_rdiag(r, len);
        A.submat(start, start-len, end, end-len) = create_rdiag(r, len);
    }
    return A;
}
arma::sp_cx_dmat Diff::create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    for (int ii = 0; ii<len; ii++) {
        temp(ii, ii) = a(pair_to_single(ii, i, len));
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }
    return temp;
}

arma::sp_cx_dmat Diff::create_rdiag(const std::complex<double> r, const int len) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    temp.diag().fill(r);
    return temp;
}


Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    M = std::round(1./h);
    len = (M-2);
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    this->u = wave_init(xc, yc, widthx, widthy, px, py, M);
    AB = create_AB(V, h, dt, M);
    n_timesteps = std::round(T/dt);
    storage = arma::cube(len, len, n_timesteps + 1);
    Re = arma::cube(len, len, n_timesteps + 1);
    Im = arma::cube(len, len, n_timesteps + 1);
    this->T = T;
    this->dt = dt;
}

void Diff::run() {
    storage.slice(0) = this->probability(this->u);
    print_u(0);
    for (int i=1; i < n_timesteps + 1 ; i++) {
        b = AB.at(1) * this->u;
        this->u = arma::spsolve(AB.at(0), b);
        storage.slice(i) = this->probability(this->u);
        print_u(i);
    }
}



arma::mat Diff::probability(arma::cx_dvec &u) {
    arma::mat prob = arma::mat(len, len);
    for (int j=0; j < len; j++) {
        for (int i=0; i < len; i++) {
            prob(i, j) = std::real(std::conj(this->u(Diff::pair_to_single(i, j, len))) * this->u(Diff::pair_to_single(i, j, len)));
        }
    }
    return prob;
}

void Diff::print_u(int t) {

    arma::cx_mat temp = arma::cx_mat(len,len);
    int temp_len = len*len;
    for (int i = 0; i < temp_len; i++) {
        std::tuple<int,int> a = single_to_pair(i,len);
        temp(std::get<0>(a),std::get<1>(a)) = this->u(i);
        }
    Re.slice(t) = arma::real(temp);
    Im.slice(t) = arma::imag(temp);

}


arma::cx_dvec Diff::wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
    int len = M-2;
    arma::cx_dvec u = arma::cx_dvec(len*len);
    double h = 1./len;
    std::complex<double> im(0., 1.);
    for (int i=0; i<(len); i++) {
        for (int j=0; j<(len); j++) {
            double temp_x = ((h*(j+1)) - centerx);
            double temp_y = ((h*(i+1)) - centery);
            std::complex<double> arg = std::exp(-(temp_x*temp_x/(2*widthx*widthx)) - (temp_y*temp_y/(2*widthy*widthy)) + im*(px*temp_x) + im*(py*temp_y));
            u(pair_to_single(i, j, len)) = arg;
        }
    }
    u = u / std::sqrt(arma::accu(arma::conj(u) % u));

    return u;
}


arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M-2;
    double h = 1./len;
    arma::mat V = arma::mat(len, len).zeros();
    int j_start = std::round((centerx - (thickx/2))/h);
    int j_end = std::round((centerx + (thickx/2))/h) - 1;
    bool center_open = n_slit%2 != 0;   

    int counter = n_slit;
    int i_center = std::round(0.5/h) -1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture/h);
    int len_slit_sep = std::round(slit_sep/h);

    if (center_open) {
        int i_loc_up = i_center - std::round(aperture/(2*h));
        int i_loc_down = i_center + std::round(aperture/(2*h));
        i_loc_up--;
        i_loc_down++;
        counter--; 
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }

    else if (n_slit == 0) {
        for (int j=j_start; j<j_end; j++) {
            V.col(j).fill(potential);
        }
        return V;
    }
    else {
        int temp = std::round(len_slit_sep/2.);
        for (int ii=0; ii<temp; ii++) {
            for (int j=j_start; j<j_end; j++) {
                V(i_loc_up, j) = potential;
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }

        return V;
    }
}*/

/*std::vector<arma::sp_cx_dmat> Diff::create_AB(arma::mat &V, const double h, const double dt, const int M) {
    int len = M - 2;
    std::complex<double> im(0., 1.);
    std::complex<double> r = im * dt / (2 * h * h);
    arma::cx_dvec a = arma::cx_dvec(len * len);
    arma::cx_dvec b = arma::cx_dvec(len * len);
    for (int ii = 0; ii < len; ii++) {
        arma::vec temp_col = V.col(ii);
        for (int i = 0; i < len; i++) {
            std::complex<double> temp_im = im * ((dt * temp_col(i)) / 2.);
            a(ii * len + i) = 1. + 4. * r + temp_im;
            b(ii * len + i) = 1. - 4. * r - temp_im;
        }
    }
    return std::vector<arma::sp_cx_dmat>{create_mat(a, -r, len), create_mat(b, r, len)};
}

arma::sp_cx_dmat Diff::create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len) {
    int lenlen = len * len;
    arma::sp_cx_dmat A = arma::sp_cx_dmat(lenlen, lenlen);
    for (int i = 0; i < len; i++) {
        int start = len * i;
        int end = len * i + (len - 1);
        A.submat(start, start, end, end) = create_tri(a, r, len, i);
    }
    for (int i = 1; i < len; i++) {
        int start = len * i;
        int end = len * i + (len - 1);
        A.submat(start - len, start, end - len, end) = create_rdiag(r, len);
        A.submat(start, start - len, end, end - len) = create_rdiag(r, len);
    }
    return A;
}
arma::sp_cx_dmat Diff::create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    for (int ii = 0; ii < len; ii++) {
        // Directly compute 1D index
        int index = ii + i * len;
        temp(ii, ii) = a(index);
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }
    return temp;
}

arma::sp_cx_dmat Diff::create_rdiag(const std::complex<double> r, const int len) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    temp.diag().fill(r);
    return temp;
}


Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    M = std::round(1./h);
    len = (M-2);
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    this->u = wave_init(xc, yc, widthx, widthy, px, py, M);
    AB = create_AB(V, h, dt, M);
    n_timesteps = std::round(T/dt);
    storage = arma::cube(len, len, n_timesteps + 1);
    Re = arma::cube(len, len, n_timesteps + 1);
    Im = arma::cube(len, len, n_timesteps + 1);
    this->T = T;
    this->dt = dt;
}

void Diff::run() {
    storage.slice(0) = this->probability(this->u);
    print_u(0);
    for (int i=1; i < n_timesteps + 1 ; i++) {
        b = AB.at(1) * this->u;
        this->u = arma::spsolve(AB.at(0), b);
        storage.slice(i) = this->probability(this->u);
        print_u(i);
    }
}



arma::mat Diff::probability(arma::cx_dvec &u) {
    arma::mat prob = arma::mat(len, len);
    for (int k = 0; k < len * len; k++) {
        // Directly compute 1D indices
        int i = k % len;
        int j = k / len;
        prob(i, j) = std::real(std::conj(u(k)) * u(k));
    }
    return prob;
}

void Diff::print_u(int t) {
    arma::cx_mat temp = arma::cx_mat(len, len);
    int temp_len = len * len;
    for (int i = 0; i < temp_len; i++) {
        // Directly compute 2D indices
        int row = i % len;
        int col = i / len;
        temp(row, col) = u(i);
    }
    Re.slice(t) = arma::real(temp);
    Im.slice(t) = arma::imag(temp);
}


arma::cx_dvec Diff::wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
    int len = M - 2;
    arma::cx_dvec u = arma::cx_dvec(len * len);
    double h = 1.0 / len;
    std::complex<double> im(0., 1.);
    for (int i = 0; i < len; i++) {
        for (int j = 0; j < len; j++) {
            double temp_x = ((h * (j + 1)) - centerx);
            double temp_y = ((h * (i + 1)) - centery);
            std::complex<double> arg = std::exp(-(temp_x * temp_x / (2 * widthx * widthx)) - (temp_y * temp_y / (2 * widthy * widthy)) + im * (px * temp_x) + im * (py * temp_y));
            u(i * len + j) = arg;
        }
    }
    u = u / std::sqrt(arma::accu(arma::conj(u) % u));

    return u;
}

arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M - 2;
    double h = 1.0 / len;
    arma::mat V = arma::mat(len, len).zeros();
    int j_start = std::round((centerx - (thickx / 2)) / h);
    int j_end = std::round((centerx + (thickx / 2)) / h) - 1;
    bool center_open = n_slit % 2 != 0;

    int counter = n_slit;
    int i_center = std::round(0.5 / h) - 1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture / h);
    int len_slit_sep = std::round(slit_sep / h);

    if (center_open) {
        int i_loc_up = i_center - std::round(aperture / (2 * h));
        int i_loc_down = i_center + std::round(aperture / (2 * h));
        i_loc_up--;
        i_loc_down++;
        counter--;
        while (counter > 0) {
            for (int ii = 0; ii < len_slit_sep; ii++) {
                for (int j = j_start; j < j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i = i_loc_up; i >= 0; i--) {
            for (int j = j_start; j < j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i = i_loc_down; i < len; i++) {
            for (int j = j_start; j < j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    } else if (n_slit == 0) {
        for (int j = j_start; j < j_end; j++) {
            V.col(j).fill(potential);
        }
        return V;
    } else {
        int temp = std::round(len_slit_sep / 2.0);
        for (int ii = 0; ii < temp; ii++) {
            for (int j = j_start; j < j_end; j++) {
                V(i_loc_up, j) = potential;
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;
        while (counter > 0) {
            for (int ii = 0; ii < len_slit_sep; ii++) {
                for (int j = j_start; j < j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i = i_loc_up; i >= 0; i--) {
            for (int j = j_start; j < j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i = i_loc_down; i < len; i++) {
            for (int j = j_start; j < j_end; j++) {
                V(i, j) = potential;
            }
        }

        return V;
    }
}*/


/*int Diff::pair_to_single(const int i, const int j, const int len) {
    return (j) * len + (i);
}

std::tuple<int, int> Diff::single_to_pair(const int k, const int len) {
    int i = k % len;
    int j = k / len;
    return std::tuple<int, int>{i, j};
}

std::vector<arma::sp_cx_dmat> Diff::create_AB(arma::mat &V, const double h, const double dt, const int M) {
    int len = M-2;
    std::complex<double> im(0., 1.);
    std::complex<double> r = im*dt/(2*h*h);
    arma::cx_dvec a = arma::cx_dvec(len*len);
    arma::cx_dvec b = arma::cx_dvec(len*len);
    for (int ii=0; ii < (len); ii++) {
        arma::vec temp_col = V.col(ii);
        for (int i=0; i < (len); i++) {
            std::complex<double> temp_im = im * ((dt*temp_col(i)) / 2.);
            a(pair_to_single(i, ii, len)) = 1. + 4.*r + temp_im;
            b(pair_to_single(i, ii, len)) = 1. - 4.*r - temp_im;
        }
    }
    return std::vector<arma::sp_cx_dmat>{create_mat(a, -r, len), create_mat(b, r, len)};
}

arma::sp_cx_dmat Diff::create_mat(arma::cx_dvec &a, const std::complex<double> r, const int len) {
    int lenlen = len*len;
    arma::sp_cx_dmat A = arma::sp_cx_dmat(lenlen, lenlen);
    for (int i = 0; i < len; i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start, start, end, end) = create_tri(a, r, len, i);
    }
    for (int i = 1; i<(len); i++) {
        int start = len*i;
        int end = len*i + (len-1);
        A.submat(start-len, start, end-len, end) = create_rdiag(r, len);
        A.submat(start, start-len, end, end-len) = create_rdiag(r, len);
    }
    return A;
}
arma::sp_cx_dmat Diff::create_tri(arma::cx_dvec &a, const std::complex<double> r, const int len, const int i) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    for (int ii = 0; ii<len; ii++) {
        temp(ii, ii) = a(pair_to_single(ii, i, len));
        temp.diag(1).fill(r);
        temp.diag(-1).fill(r);
    }
    return temp;
}

arma::sp_cx_dmat Diff::create_rdiag(const std::complex<double> r, const int len) {
    arma::sp_cx_dmat temp = arma::sp_cx_dmat(len, len);
    temp.diag().fill(r);
    return temp;
}


Diff::Diff(double h, double dt, double T, double xc, double yc, double px, double py, double widthx, double widthy, double potential, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    M = std::round(1./h);
    len = (M-2);
    V = this->potential(potential, M, n_slit, thickx, centerx, slit_sep, aperture);
    this->u = wave_init(xc, yc, widthx, widthy, px, py, M);
    AB = create_AB(V, h, dt, M);
    n_timesteps = std::round(T/dt);
    storage = arma::cube(len, len, n_timesteps + 1);
    Re = arma::cube(len, len, n_timesteps + 1);
    Im = arma::cube(len, len, n_timesteps + 1);
    this->T = T;
    this->dt = dt;
}

void Diff::run() {
    storage.slice(0) = this->probability(this->u);
    print_u(0);
    for (int i=1; i < n_timesteps + 1 ; i++) {
        b = AB.at(1) * this->u;
        this->u = arma::spsolve(AB.at(0), b);
        storage.slice(i) = this->probability(this->u);
        print_u(i);
    }
}

void Diff::print(std::string filename) {
    storage.save(filename + ".bin");
}

void Diff::print_potential(std::string filename) {
    V.save(filename + "_pot.bin");
}

void Diff::save_u(std::string filename) {
    Re.save(filename + "_Re.bin");
    Im.save(filename + "_Im.bin");
}

arma::mat Diff::probability(arma::cx_dvec &u) {
    arma::mat prob = arma::mat(len, len);
    for (int j=0; j < len; j++) {
        for (int i=0; i < len; i++) {
            prob(i, j) = std::real(std::conj(this->u(Diff::pair_to_single(i, j, len))) * this->u(Diff::pair_to_single(i, j, len)));
        }
    }
    return prob;
}

void Diff::print_u(int t) {

    arma::cx_mat temp = arma::cx_mat(len,len);
    int temp_len = len*len;
    for (int i = 0; i < temp_len; i++) {
        std::tuple<int,int> a = single_to_pair(i,len);
        temp(std::get<0>(a),std::get<1>(a)) = this->u(i);
        }
    Re.slice(t) = arma::real(temp);
    Im.slice(t) = arma::imag(temp);

}


arma::cx_dvec Diff::wave_init(double centerx, double centery, double widthx, double widthy, double px, double py, int M) {
    int len = M-2;
    arma::cx_dvec u = arma::cx_dvec(len*len);
    double h = 1./len;
    std::complex<double> im(0., 1.);
    for (int i=0; i<(len); i++) {
        for (int j=0; j<(len); j++) {
            double temp_x = ((h*(j+1)) - centerx);
            double temp_y = ((h*(i+1)) - centery);
            std::complex<double> arg = std::exp(-(temp_x*temp_x/(2*widthx*widthx)) - (temp_y*temp_y/(2*widthy*widthy)) + im*(px*temp_x) + im*(py*temp_y));
            u(pair_to_single(i, j, len)) = arg;
        }
    }
    u = u / std::sqrt(arma::accu(arma::conj(u) % u));

    return u;
}


arma::mat Diff::potential(double potential, int M, int n_slit, double thickx, double centerx, double slit_sep, double aperture) {
    int len = M-2;
    double h = 1./len;
    arma::mat V = arma::mat(len, len).zeros();
    int j_start = std::round((centerx - (thickx/2))/h);
    int j_end = std::round((centerx + (thickx/2))/h) - 1;
    bool center_open = n_slit%2 != 0;   

    int counter = n_slit;
    int i_center = std::round(0.5/h) -1;
    int i_loc_up = i_center;
    int i_loc_down = i_center;
    int len_aperture = std::round(aperture/h);
    int len_slit_sep = std::round(slit_sep/h);

    if (center_open) {
        int i_loc_up = i_center - std::round(aperture/(2*h));
        int i_loc_down = i_center + std::round(aperture/(2*h));
        i_loc_up--;
        i_loc_down++;
        counter--; 
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        return V;
    }

    else if (n_slit == 0) {
        for (int j=j_start; j<j_end; j++) {
            V.col(j).fill(potential);
        }
        return V;
    }
    else {
        int temp = std::round(len_slit_sep/2.);
        for (int ii=0; ii<temp; ii++) {
            for (int j=j_start; j<j_end; j++) {
                V(i_loc_up, j) = potential;
                V(i_loc_down, j) = potential;
            }
            i_loc_up--;
            i_loc_down++;
        }
        i_loc_up -= len_aperture;
        i_loc_down += len_aperture;
        counter -= 2;
        while (counter > 0) {
            for (int ii=0; ii<len_slit_sep; ii++) {
                for (int j=j_start; j<j_end; j++) {
                    V(i_loc_up, j) = potential;
                    V(i_loc_down, j) = potential;
                }
                i_loc_up--;
                i_loc_down++;
            }
            i_loc_up -= len_aperture;
            i_loc_down += len_aperture;
            counter -= 2;
        }
        for (int i=i_loc_up; i>=0; i--) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }
        for (int i=i_loc_down; i<len; i++) {
            for (int j=j_start; j<j_end; j++) {
                V(i, j) = potential;
            }
        }

        return V;
    }
}*/