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

vector<string> bfs(const ChessState& start) {
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


    format_res("BFS", start.size, nodes, solutions, t1, t2);
    return res;
}