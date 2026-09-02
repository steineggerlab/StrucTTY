#pragma once
#include <unistd.h>

namespace Terminal {

struct Size {
    int rows;
    int cols;
};

struct MouseEvent {
    int  x;
    int  y;
    bool pressed;
    bool moved;
};

constexpr int KEY_UP    = 0x100;
constexpr int KEY_DOWN  = 0x101;
constexpr int KEY_LEFT  = 0x102;
constexpr int KEY_RIGHT = 0x103;
constexpr int KEY_MOUSE = 0x104;

void enter_raw_mode();
void exit_raw_mode();

Size get_size();

void clear();
void cursor_home();
void move_cursor(int row, int col);
void clear_to_eol();

void flush_input();

int  read_key();

bool read_mouse(MouseEvent& ev);

}
