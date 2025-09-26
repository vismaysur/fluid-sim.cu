#include "raylib.h"
#include "graphics.h"
#include <cstring>

#define N 20
#define IDX(N_, i, j) (i * N_ + j)
#define U_IDX(N_, i, j) (i * (N_ + 1) + j)
#define V_IDX(N_, i, j) (i * N_ + j)

float rho[(N + 2) * (N + 2)];
float u[(N + 3) * (N + 2)];
float v[(N + 2) * (N + 3)];

float sources[(N + 2) * (N + 2)];

int initialize(int N_, float *rho_)
{
    int size = (N_ + 2) * (N_ + 2); // + 2 for ghost cells / boundary

    for (int i = 0; i < size; i++)
    {
        rho_[i] = 0.0f;
    }

    return 0;
}

int addSource(int N_, float *rho_, float *source, float dt)
{
    int size = (N_ + 2) * (N_ + 2); // + 2 for ghost cells / boundary

    for (int i = 0; i < size; i++)
    {
        rho_[i] += dt * source[i];
    }

    return 0;
}

int main()
{
    InitWindow(WINDOW_DIM, WINDOW_DIM, "Fluid-Simulation");
    SetTargetFPS(60);

    initialize(N, rho);

    while (!WindowShouldClose())
    {
        memset(sources, 0, (N + 2) * (N + 2) * sizeof(float));

        handleMouseInput(sources, N);

        addSource(N, rho, sources, 0.016f);

        // diffuse

        // advect

        drawGrid(rho, N);
    }

    CloseWindow();
}