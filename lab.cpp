#include <bits/stdc++.h>
using namespace std;

double f(double x) { return x*x + 5*sin(x) - 1; }
double df(double x) { return 2*x + 5*cos(x); }

pair<bool,double> phi_and_ok(double x) {
    if (x == 0) return {false, 0.0};
    double val = x + (x*x + 5*sin(x) - 1) / (15*x*x);
    return {true, val};
}

pair<bool,double> phi_prime_and_ok(double x) {
    if (x == 0) return {false, 0.0};
    double val = 1 + (-10*sin(x) + 5*x*cos(x) + 2) / (15*x*x*x);
    return {true, val};
}

vector<double> linspace(double a, double b, int n) {
    vector<double> v; v.reserve(n);
    for (int i=0;i<n;i++){
        double t = (n==1) ? 0.0 : double(i)/(n-1);
        v.push_back(a + t*(b-a));
    }
    return v;
}

double simpleIteration(double x0, double a, double b, double q_phi, double eps, int maxIter) {
    vector<double> iter_vals;
    iter_vals.push_back(x0); 

    double q = q_phi;
    double post_est = (1 - q)/q * eps;

    int n_aprior_phi = floor(log(fabs(phi_and_ok(x0).second - x0) / ((1 - q)*eps)) / log(1/q)+1);
    int N = min(n_aprior_phi, maxIter);

    cout << "Апрiорна оцiнка кiлькостi iтерацiй: n >= " << n_aprior_phi << "\n\n";
    cout << setw(4) << "k" << setw(20) << "x_k" << setw(20) << "|x_k - x_{k-1}|\n";
    cout << string(45,'-') << "\n";
    cout << setw(4) << 0 << setw(20) << iter_vals[0] << setw(20) << "-" << "\n";

    double x_prev = x0, x_curr;
    for(int k=1; k<=N; k++){
        x_curr = phi_and_ok(x_prev).second;
        iter_vals.push_back(x_curr);
        cout << setw(4) << k << setw(20) << x_curr << setw(20) << fabs(x_curr - x_prev) << "\n";
        x_prev = x_curr;
    }

    double x_final = iter_vals.back();
    cout << "\n Корiнь ≈ " << x_final << "\n";

    return x_final;
}

double relaxationMethod(double x0, double a, double b, const vector<double>& pts, double eps, int maxIter) {
    double m1 = numeric_limits<double>::infinity(), M1 = 0.0;
    for(double x: pts){
        double val = fabs(df(x));
        m1 = min(m1,val);
        M1 = max(M1,val);
    }
    double tau = 2.0/(M1 + m1);
    double q_relax = (M1 - m1)/(M1 + m1);

    double delta = max(fabs(x0-a), fabs(x0-b));
    int n_aprior_relax = max(1,(int)ceil(log(delta/eps)/log(1.0/q_relax)));
    
    cout << "m1 = " << m1 << ", M1 = " << M1 << "\n";
    cout << "Оптимальний tau ≈ " << tau << ", q_relax ≈ " << q_relax << "\n";
    cout << "Апрiорна оцiнка кiлькостi iтерацiй: n >= " << n_aprior_relax << "\n\n";

    cout << setw(4) << "k" << setw(20) << "x_k" << setw(20) << "|x_k - x_{k-1}|\n";
    cout << string(45,'-') << "\n";

    double xr_prev = x0, xr_curr;
    cout << setw(4) << 0 << setw(20) << xr_prev << setw(20) << "-" << "\n";

    for(int k=1; k<=n_aprior_relax; k++){
        xr_curr = xr_prev + tau*f(xr_prev);
        cout << setw(4) << k << setw(20) << xr_curr << setw(20) << fabs(xr_curr - xr_prev) << "\n";
        xr_prev = xr_curr;
    }

    cout << "\nКорiнь ≈ " << xr_prev << "\n";

    return xr_prev;
}


int main() {
    ios::sync_with_stdio(false); cin.tie(nullptr);
    cout << fixed << setprecision(4);

    const double eps = 1e-4;
    const int maxIter = 1000;

    cout << "=== Перевiрка значень f(x) на контрольних точках ===\n";
    vector<double> checks = {0.0, -1.0, -2.0, -3.0};
    for(double x: checks) cout << "f(" << x << ") = " << f(x) << "\n";

    cout << "\nОскiльки f(-3) > 0 i f(-2) < 0, маємо змiну знака на промiжку (-3, -2).\n";
    cout << "Отже, шукаємо вiд'ємний корiнь в I0 = (-3, -2).\n\n";

    double a = -2.3, b = -2.1, x0 = -2.25;
    cout << "Обрано вузький iнтервал I = [" << a << ", " << b << "] та початкове наближення x0 = " << x0 << "\n\n";

    int samples = 400;
    auto pts = linspace(a,b,samples);
    double q_phi = 0.0;
    for(double x: pts){
        auto pr = phi_prime_and_ok(x);
        if(!pr.first){ cout << "phi не визначена на I\n"; return 0; }
        q_phi = max(q_phi, fabs(pr.second));
    }
    cout << "Перевiрка phi на I: max |phi'(x)| ≈ " << q_phi << "\n";
    if(q_phi >= 1.0){ cout << "q >=1, метод не сходиться на iнтервалi\n"; return 0; }
    cout << "Оскiльки q < 1, метод простої iтерації може сходитися.\n\n";

    cout << "=== Метод простої iтерацiї ===\n";
    double x_iter = simpleIteration(x0, a, b, q_phi, eps, maxIter);

    cout << "=== Метод релаксацiї ===\n";
    double x_relax = relaxationMethod(x0, a, b, pts, eps, maxIter);

    cout << "=== Порiвняння ===\n";
    cout << "Метод простої iтерацiї: x ≈ " << x_iter << "\n";
    cout << "Метод релаксацiї      : x ≈ " << x_relax << "\n";

    return 0;
}


