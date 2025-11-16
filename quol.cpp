#include <iostream>
#include <complex>
#include <vector>
#include <iomanip>
#include <map>
#include <numbers>

using namespace std;
const double PI = numbers::pi;

typedef complex<double> CD;

struct QFTResult {
    map<int, CD> state;

    QFTResult() {}
    QFTResult(map<int, CD> m) : state(m) {}

    QFTResult operator+(const QFTResult& rhs) const {
        map<int, CD> res = state;
        for (auto &[num, xy] : rhs.state) {
            res[num] += xy;
        }
        return QFTResult(res);
    }

    void normalize() {
        double totalNorm = 0;
        for (auto &[_, amp] : state) totalNorm += norm(amp);
        if (totalNorm > 0) {
            double factor = 1.0 / sqrt(totalNorm);
            for (auto &[_, amp] : state) amp *= factor;
        }
    }
};

struct QFT {
    int m;
    CD ampOfOne;

    QFT(int m) : m(m) {
        ampOfOne = CD(1.0 / sqrt(m), 0);
    }

    QFTResult qft(int x) const {
        QFTResult res;
        for (int y = 0; y < m; y++)
            res.state[y] = ampOfOne * polar(1.0, 2 * PI * x * y / m);
        return res;
    }

    void print(const QFTResult& res) const {
        bool first = true;
        for (auto &[num, amp] : res.state) {
            double p = norm(amp);
            if (p < 0.005) continue;
            if (!first) cout << " + ";
            cout << fixed << setprecision(4) << p << "|" << num << ">";
            first = false;
        }
    }
};

int main() {
    int n, a, m;
    cin >> n >> a >> m;

    map<int, vector<int>> M;
    int aToK = 1;
    M[aToK].push_back(0);
    for (int k = 1; k < m; k++) {
        (aToK *= a) %= n;
        M[aToK].push_back(k);
    }

    QFT qft(m);

    for (auto &[right, nums] : M) {
        QFTResult res;
        for (auto &num : nums)
            res = res + qft.qft(num);

        res.normalize();
        cout << "( ";
        qft.print(res);
        cout << " ) |" << right << ">\n";
    }
}
