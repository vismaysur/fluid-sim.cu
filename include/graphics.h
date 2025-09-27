#include "raylib.h"
#include <iostream>

#define IDX(N_, i, j) ((j) * (N_) + (i))
#define WINDOW_DIM 800.0f

// Functions assume (grid_dim * grid_dim) dimensions

void handleMouseInput(float *sources, int grid_dim)
{
    Vector2 mouse = GetMousePosition();

    int grid_x = (int)(mouse.x / WINDOW_DIM * (grid_dim - 2));
    int grid_y = (int)(mouse.y / WINDOW_DIM * (grid_dim - 2));

    if (grid_x >= 0 && grid_x < grid_dim - 1 && grid_y >= 0 && grid_y < grid_dim - 1)
    {
        int idx = IDX(grid_dim, grid_x + 1, grid_y + 1);

        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON))
        {
            sources[idx] = 5.0f;
        }

        if (IsMouseButtonDown(MOUSE_RIGHT_BUTTON))
        {
            sources[idx] = -5.0f;
        }
    }
}

void drawGrid(float *density, int grid_dim)
{
    BeginDrawing();

    ClearBackground(WHITE);

    int cell_size = WINDOW_DIM / (grid_dim - 2);

    for (int grid_x = 1; grid_x < grid_dim - 1; grid_x++)
    {
        for (int grid_y = 1; grid_y < grid_dim - 1; grid_y++)
        {
            int idx = IDX(grid_dim, grid_x, grid_y);

            float d = density[idx];

            if (d >= 0.01f)
            {
                Color color = {0, 0, 255, static_cast<unsigned char>(d * 255)};

                DrawRectangle(
                    (grid_x - 1) * cell_size,
                    (grid_y - 1) * cell_size,
                    cell_size,
                    cell_size,
                    color);
            }
        }
    }

    EndDrawing();
}