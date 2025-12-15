#include <iostream>
#include <vector>
#include <queue>
#include <unordered_set>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <random>

using namespace std;

struct ChessState {
    uint8_t board_size;
    uint8_t depth;
    vector<int8_t> positions;

    uint64_t used_cols;
    uint64_t used_main_diag;
    uint64_t used_side_diag;

    explicit ChessState(uint8_t n)
        : board_size(n),
        depth(0),
        positions(n, -1),
        used_cols(0),
        used_main_diag(0),
        used_side_diag(0) {
    }

    inline uint64_t main_diag(int r, int c) const {
        return r - c + board_size - 1;
    }

    inline uint64_t side_diag(int r, int c) const {
        return r + c;
    }

    bool can_place(int r, int c) const {
        return !(
            (used_cols & (1ULL << c)) ||
            (used_main_diag & (1ULL << main_diag(r, c))) ||
            (used_side_diag & (1ULL << side_diag(r, c)))
            );
    }

    void apply_move(int r, int c) {
        positions[r] = c;
        used_cols |= (1ULL << c);
        used_main_diag |= (1ULL << main_diag(r, c));
        used_side_diag |= (1ULL << side_diag(r, c));
        ++depth;
    }

    void rollback_move(int r, int c) {
        positions[r] = -1;
        used_cols &= ~(1ULL << c);
        used_main_diag &= ~(1ULL << main_diag(r, c));
        used_side_diag &= ~(1ULL << side_diag(r, c));
        --depth;
    }

    bool is_complete() const {
        return depth == board_size;
    }

    string to_string() const {
        string out;
        for (int i = 0; i < board_size; ++i) {
            out += (positions[i] < 0 ? '.' : char('A' + positions[i]));
        }
        return out;
    }

    bool operator==(const ChessState& other) const {
        return positions == other.positions;
    }
};

template<>
struct hash<ChessState> {
    size_t operator()(const ChessState& s) const noexcept {
        size_t h = 0;
        for (auto p : s.positions)
            h = h * 31 + (p + 1);
        return h;
    }
};

// BFS
pair<vector<string>, pair<int, int>> bfs(const ChessState& start) {
    queue<ChessState> q;
    unordered_set<ChessState> visited;
    int nodes = 0, solutions = 0;
    vector<string> res;

    q.push(start);
    visited.insert(start);

    while (!q.empty()) {
        ChessState cur = q.front(); q.pop();
        nodes++;

        if (cur.is_complete()) {
            solutions++;
            res.push_back(cur.to_string());
            continue;
        }

        for (int col = 0; col < cur.board_size; col++) {
            if (cur.can_place(cur.depth, col)) {
                ChessState next = cur;
                next.apply_move(cur.depth, col);
                if (!visited.count(next)) {
                    visited.insert(next);
                    q.push(next);
                }
            }
        }
    }

    return { res, {nodes, solutions} };
}

// DFS
pair<vector<string>, pair<int, int>> dfs(const ChessState& start) {
    int nodes = 0;
    int solutions = 0;
    vector<string> res;
    ChessState state = start;

    auto dfs_impl = [&](auto&& self, int row) -> void {
        ++nodes;

        if (state.is_complete()) {
            ++solutions;
            res.push_back(state.to_string());
            return;
        }

        for (int col = 0; col < state.board_size; ++col) {
            if (state.can_place(row, col)) {
                state.apply_move(row, col);
                self(self, row + 1);
                state.rollback_move(row, col);
            }
        }
    };

    dfs_impl(dfs_impl, 0);

    return { res, {nodes, solutions} };
}

// IDS
pair<vector<string>, pair<int, int>> ids(const ChessState& start) {
    int nodes = 0;
    int solutions = 0;
    vector<string> res;
    ChessState state = start;

    auto ids_impl = [&](auto&& self, int row, int limit) -> void {
        ++nodes;

        if (state.is_complete()) {
            ++solutions;
            res.push_back(state.to_string());
        }

        if (row != limit) {
            for (int col = 0; col < state.board_size; ++col) {
                if (state.can_place(row, col)) {
                    state.apply_move(row, col);
                    self(self, row + 1, limit);
                    state.rollback_move(row, col);
                }
            }
        }
        };

    for (int depth = 1; depth <= start.board_size; ++depth) {
        ids_impl(ids_impl, 0, depth);
    }

    return { res, {nodes, solutions} };
}

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

int rand_int(int a, int b) { // [a,b]
    uniform_int_distribution<int> dist(a, b);
    return dist(rng);
}

double rand01() {
    uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
}

int count_conflicts(const vector<int>& queens) {
    int N = queens.size();
    vector<int> col_count(N, 0);
    vector<int> diag1_count(2 * N - 1, 0);
    vector<int> diag2_count(2 * N - 1, 0);

    for (int row = 0; row < N; row++) {
        int c = queens[row];
        col_count[c]++;
        diag1_count[row - c + N - 1]++;
        diag2_count[row + c]++;
    }

    int conflicts = 0;
    auto add_conf = [&](const vector<int>& v) { // для подсчета количества конфликтных пар
        for (int cnt : v)
            if (cnt > 1) conflicts += cnt * (cnt - 1) / 2;
        };

    add_conf(col_count);
    add_conf(diag1_count);
    add_conf(diag2_count);

    return conflicts;
}

vector<int> ga_nqueens(int N, int pop_size, int generations, double mutation_rate) {

    // На этом этапе формируется начальная популяция генетического алгоритма. 
    // Каждая особь представляется перестановкой столбцов, 
    // что гарантирует отсутствие конфликтов по строкам и столбцам.
    vector<vector<int>> population;
    for (int i = 0; i < pop_size; i++) {
        vector<int> q(N);
        for (int j = 0; j < N; j++) q[j] = j;
        shuffle(q.begin(), q.end(), rng);
        population.push_back(q);
    }

    // Один проход = одно поколение генетического алгоритма.
    for (int gen = 0; gen < generations; gen++) {
        
        // (количество конфликтов, индекс особи)
        vector<pair<int, int>> fitness;
        
        for (int i = 0; i < pop_size; i++)
            fitness.push_back({ count_conflicts(population[i]), i }); // Меньше конфликтов -> лучше особь
        sort(fitness.begin(), fitness.end()); // Сортировка по возрастанию конфликтов

        if (fitness[0].first == 0) return population[fitness[0].second];

        // Новая популяция
        vector<vector<int>> new_pop;
        for (int k = 0; k < pop_size; k++) {
            // Выбор родителей, рандомно из второй половины популяции
            int i1 = fitness[rand_int(0, pop_size / 2)].second;
            int i2 = fitness[rand_int(0, pop_size / 2)].second;

            vector<int> child(N, -1);
            int cut = rand_int(0, N - 1); // разделение особи

            for (int j = 0; j <= cut; j++) child[j] = population[i1][j]; // Копирование части от первого родителя
            
            // столбцы, уже попавшие в потомка
            vector<bool> used(N, false);
            for (int j = 0; j <= cut; j++) used[child[j]] = true;

            // Добавление генов второго родителя
            int idx = cut + 1;
            for (int j = 0; j < N; j++) {
                if (!used[population[i2][j]]) child[idx++] = population[i2][j];
            }

            // С вероятностью mutation_rate выполняется мутация,
            // случайно меняются два гена местами
            if (rand01() < mutation_rate) {
                int a = rand_int(0, N - 1), b = rand_int(0, N - 1);
                swap(child[a], child[b]); 
            }

            new_pop.push_back(child);
        }
        // Прошлое поколение заменяется потомками 
        population = new_pop;
    }

    return {};
}

string genetic_algorithm(int N, int pop_size = 100, int generations = 1000, double mutation_rate = 0.1) {

    vector<int> solution;

    while (solution.empty())
        solution = ga_nqueens(N, pop_size, generations, mutation_rate);

    string result;

    if (!solution.empty()) {
        for (int col : solution)
            result += char('A' + col);
    }

    return result;
}

// Simulated Annealing
vector<int> simulated_annealing_nqueens(int N, double t_start = 1000.0, double t_end = 1e-4, double cooling_rate = 0.995) {
    // Инициализация начального состояния
    vector<int> current(N);
    for (int i = 0; i < N; i++) current[i] = i;
    shuffle(current.begin(), current.end(), rng);

    // Оценка состояния
    int current_conflicts = count_conflicts(current);
    vector<int> best = current; 
    int best_conflicts = current_conflicts;

    double t = t_start;

    while (t > t_end && best_conflicts > 0) {
        // Меняем местами два случайных ферзя, те выбираем случайного "соседа"
        vector<int> next = current;
        int a = rand_int(0, N - 1);
        int b = rand_int(0, N - 1);
        swap(next[a], next[b]);

        int next_conflicts = count_conflicts(next);

        // разница между новым и текущим состоянием
        // delta_e < 0 -> новое состояние лучше;
        // delta_e > 0 -> новое состояние хуже.
        int delta_e = next_conflicts - current_conflicts;
        
        // Если новое состояние лучше (delta_e < 0) -> всегда принимаем.
        // Если хуже (delta_e > 0) -> принимаем с вероятностью exp(-delta_e / t).
        // Вероятность уменьшается с температурой.
        // Позволяет выпрыгнуть из состояния локального минимума
        if (delta_e < 0 || rand01() < exp(-delta_e / t)) {
            current = next;
            current_conflicts = next_conflicts;

            if (current_conflicts < best_conflicts) {
                best = current;
                best_conflicts = current_conflicts;
            }
        }
        t *= cooling_rate;
    }

    if (best_conflicts == 0) return best;
    else return {};
}

string simulated_annealing(int N) {
    vector<int> solution;

    while (solution.empty())
        solution = simulated_annealing_nqueens(N);

    string result;
    for (int col : solution)
        result += char('A' + col);

    return result;
}

int main() {
    setlocale(LC_ALL, "Russian");

    int N;
    cout << "Введите размер доски N: ";
    cin >> N;

    ChessState start(N);

    auto bfs_res = ids(start);

    for (auto item : bfs_res.first) {
        cout << item << endl;
    }

    cout << bfs_res.second.second << endl;
}