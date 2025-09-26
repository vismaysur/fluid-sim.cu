#include "raylib.h"
#include <iostream>

#define WINDOW_DIM 800.0f

// Functions assume (grid_dim * grid_dim) dimensions

void handleMouseInput(float *sources, int grid_dim)
{
    Vector2 mouse = GetMousePosition();

    int grid_x = (int)(mouse.x / WINDOW_DIM * grid_dim);
    int grid_y = (int)(mouse.y / WINDOW_DIM * grid_dim);
    grid_y = grid_dim - grid_y - 1; // flip orientation for simulation

    if (grid_x >= 0 && grid_x < grid_dim && grid_y >= 0 && grid_y < grid_dim)
    {
        int idx = (grid_y + 1) * (grid_dim + 2) + (grid_x + 1); // Account for ghost / boundary cells.

        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON))
        {
            sources[idx] = 500.0f;
        }

        if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON))
        {
            sources[idx] = -500.0f;
        }
    }
}

void drawGrid(float *density, int grid_dim)
{
    BeginDrawing();

    ClearBackground(WHITE);

    float cell_size = WINDOW_DIM / grid_dim;

    for (int i = 1; i <= grid_dim; i++)
    {
        for (int j = 1; j <= grid_dim; j++)
        {
            int idx = i * (grid_dim + 2) + j;

            float d = density[idx] / 100.0f;

            if (d >= 0.01f)
            {
                float intensity = d > 1.0 ? 1.0 : d;

                Color color = {0, 0, 255, static_cast<unsigned char>(intensity * 255)};

                DrawRectangle(
                    (j - 1) * cell_size,
                    (grid_dim - i) * cell_size,
                    cell_size,
                    cell_size,
                    color);
            }
        }
    }

    EndDrawing();
}