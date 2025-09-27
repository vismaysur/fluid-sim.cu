#include "raylib.h"
#include "graphics.h"
#include <cstring>

#define N 40
#define IDX(N_, i, j) ((j) * (N_) + (i))
#define U_IDX(N_, i, j) ((j) * ((N_) + 1) + (i))
#define V_IDX(N_, i, j) ((j) * (N_) + (i))

float rho[N * N];
float rho_prev[N * N];
float u[(N + 1) * N];
float u_prev[(N + 1) * N];
float v[N * (N + 1)];
float v_prev[N * (N + 1)];

float sources[N * N];

void spawnBlob(int N_, float *rho, int x1, int y1, int x2, int y2)
{
    for (int i = x1; i <= x2; i++)
    {
        for (int j = y1; j <= y2; j++)
        {
            rho[IDX(N_, i, j)] = 1.0f;
        }
    }
}

void setBnd(int N_, int b, float *A)
{
    for (int i = 1; i < N_; i++)
    {
        A[IDX(N_, i, 0)] = b == 2 ? -A[IDX(N_, i, 1)] : A[IDX(N_, i, 1)];
        A[IDX(N_, i, N_ - 1)] = b == 2 ? -A[IDX(N_, i, N_ - 2)] : A[IDX(N_, i, N_ - 2)];
        A[IDX(N_, 0, i)] = b == 1 ? -A[IDX(N_, 1, i)] : A[IDX(N_, 1, i)];
        A[IDX(N_, N_ - 1, i)] = b == 1 ? -A[IDX(N_, N_ - 2, i)] : A[IDX(N_, N_ - 2, i)];
    }
}

void addSource(int N_, float *rho_, float *source, float dt)
{
    int size = N_ * N_;
    for (int idx = 0; idx < size; idx++)
    {
        if (source[idx] != 0.0f)
        {
            rho_[idx] = std::min(1.0f, rho_[idx] + dt * source[idx]);
        }
    }
}

void addForce(int N_, float *u, float *v, float *rho, float dt)
{
    float gravity = 9.81f * 0.008f;

    // Apply gravitation acceleration to vertical component of velocity.

    for (int i = 0; i < N_; i++)
    {
        for (int j = 0; j <= N_; j++)
        {
            v[V_IDX(N_, i, j)] += dt * gravity;
        }
    }

    setBnd(N_, 2, v);
}

// Gauss-Seidel solver for diffusion
// TODO: use Jacobi instead (better parallelization properties)
void diffuse(int N_, int b, float *curr, float *prev, float diff, float dt)
{
    float a = dt * N_ * N_ * diff;

    for (int k = 0; k < 20; k++)
    {
        for (int i = 1; i < N_ - 1; i++)
        {
            for (int j = 1; j < N_ - 1; j++)
            {
                float next = (prev[IDX(N_, i, j)] + a * (curr[IDX(N_, i + 1, j)] + curr[IDX(N_, i - 1, j)] + curr[IDX(N_, i, j + 1)] + curr[IDX(N_, i, j - 1)])) / (1 + 4 * a);
                curr[IDX(N_, i, j)] = std::min(next, 1.0f);
            }
        }

        setBnd(N_, b, curr);
    }
}

// Linear backtracing + interpolation solver for advection
void advect(int N_, int b, float *curr, float *prev, float *u, float *v, float dt)
{
    float dt0 = dt * N_;

    for (int i = 1; i < N_ - 1; i++)
    {
        for (int j = 1; j < N_ - 1; j++)
        {
            float u_interp = 0.5 * (u[U_IDX(N_, i, j)] + u[U_IDX(N_, i + 1, j)]);
            float v_interp = 0.5 * (v[V_IDX(N_, i, j)] + v[V_IDX(N_, i, j + 1)]);

            if (j == N_ - 3 && i == 20)
            {
                std::cout << "==" << std::endl;
                std::cout << v_interp << std::endl;
            }

            if (j == N_ - 2 && i == 20)
            {
                std::cout << "====" << std::endl;
                std::cout << v_interp << std::endl;
            }

            float x = i - dt0 * u_interp;
            float y = j - dt0 * v_interp;

            if (x < 1)
                x = 1;

            if (x > N_ - 2)
                x = N_ - 2;

            if (y < 1)
                y = 1;

            if (y > N_ - 2)
                y = N_ - 2;

            int i0 = (int)x;
            int i1 = i0 + 1;
            int j0 = (int)y;
            int j1 = j0 + 1;

            float s1 = x - i0;
            float s0 = 1 - s1;
            float t1 = y - j0;
            float t0 = 1 - t1;

            curr[IDX(N_, i, j)] = s0 * (t0 * prev[IDX(N_, i0, j0)] + t1 * prev[IDX(N_, i0, j1)]) +
                                  s1 * (t0 * prev[IDX(N_, i1, j0)] + t1 * prev[IDX(N_, i1, j1)]);
        }
    }

    setBnd(N_, b, curr);
}

// Project function to enforce incompressibility

int main()
{
    InitWindow(WINDOW_DIM, WINDOW_DIM, "Fluid-Simulation");
    SetTargetFPS(60);

    memset(rho, 0, N * N);

    spawnBlob(N, rho, 15, 15, 25, 25);

    float dt = 0.016f;

    // WATER
    float diff = 0.0f;
    float visc = 0.000001f;

    while (!WindowShouldClose())
    {
        memset(sources, 0, N * N * sizeof(float));

        handleMouseInput(sources, N);

        addSource(N, rho, sources, dt);

        addForce(N, u, v, rho, dt);

        // VELOCITY STEP:

        // diffuse
        // memcpy(u_prev, u, sizeof(float) * (N + 1) * N);
        // diffuse(N, 1, u, u_prev, visc, dt);

        // memcpy(v_prev, v, sizeof(float) * (N + 1) * N);
        // diffuse(N, 2, v, v_prev, visc, dt);

        // advect
        // memcpy(u_prev, u, sizeof(float) * (N + 1) * N);
        // memcpy(v_prev, v, sizeof(float) * (N + 1) * N);

        // advect(N, 1, u, u_prev, u_prev, v_prev, dt);
        // advect(N, 2, v, v_prev, u_prev, v_prev, dt);

        // DENSITY STEP:

        // diffuse
        // memcpy(rho_prev, rho, sizeof(float) * N * N);
        // diffuse(N, 0, rho, rho_prev, diff, dt);

        // advect
        memcpy(rho_prev, rho, sizeof(float) * N * N);
        advect(N, 0, rho, rho_prev, u, v, dt);

        drawGrid(rho, N);
    }

    CloseWindow();
}

// TODOS:
// - Fix diffusion so the rate of diff in each cell should depend on actual material
// - For example, you can't just diffuse velocity freely (it should depend on viscosity).
// - Make advection mass conserving
// - Add pressure projection
// - Use real units / values for coefficients and make it clear
//   how grid configuration must be adjusted accordingly