#include "Terminal.hpp"
#include <termios.h>
#include <sys/ioctl.h>
#include <sys/select.h>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <csignal>
#include <cerrno>

namespace Terminal {

static struct termios g_saved_termios;
static bool           g_raw_active = false;

static int  g_mouse_x       = 0;
static int  g_mouse_y       = 0;
static bool g_mouse_pressed = false;
static bool g_mouse_moved   = false;

static void do_exit_raw() {
    if (!g_raw_active) return;
    g_raw_active = false;

    static const char disable[] =
        "\033[?1049l"
        "\033[?1003l"
        "\033[?1002l"
        "\033[?1000l"
        "\033[?1006l"
        "\033[?25h";
    write(STDOUT_FILENO, disable, sizeof(disable) - 1);

    tcsetattr(STDIN_FILENO, TCSAFLUSH, &g_saved_termios);
}

static void atexit_handler()  { do_exit_raw(); }
static void sig_handler(int)  { do_exit_raw(); _exit(0); }

void enter_raw_mode() {
    if (g_raw_active) return;

    if (tcgetattr(STDIN_FILENO, &g_saved_termios) == -1) return;
    g_raw_active = true;

    struct termios raw = g_saved_termios;
    raw.c_iflag &= ~(unsigned)(BRKINT | ICRNL | INPCK | ISTRIP | IXON);
    raw.c_oflag &= ~(unsigned)(OPOST);
    raw.c_cflag |=  (unsigned)(CS8);
    raw.c_lflag &= ~(unsigned)(ECHO | ICANON | IEXTEN | ISIG);
    raw.c_cc[VMIN]  = 1;
    raw.c_cc[VTIME] = 0;
    tcsetattr(STDIN_FILENO, TCSAFLUSH, &raw);

    static const char enable[] =
        "\033[?1000h"
        "\033[?1002h"
        "\033[?1003h"
        "\033[?1006h"
        "\033[?1049h"
        "\033[?25l";
    write(STDOUT_FILENO, enable, sizeof(enable) - 1);
    fflush(stdout);

    atexit(atexit_handler);
    signal(SIGTERM, sig_handler);
    signal(SIGINT,  sig_handler);
}

void exit_raw_mode() {
    do_exit_raw();
}

Size get_size() {
    struct winsize ws;
    if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == -1 || ws.ws_row == 0) {
        return {24, 80};
    }
    return {(int)ws.ws_row, (int)ws.ws_col};
}

void clear() {
    write(STDOUT_FILENO, "\033[2J\033[H", 7);
    fflush(stdout);
}

void cursor_home() {
    write(STDOUT_FILENO, "\033[H", 3);
}

void move_cursor(int row, int col) {
    char buf[24];
    int n = snprintf(buf, sizeof(buf), "\033[%d;%dH", row + 1, col + 1);
    if (n > 0) fwrite(buf, 1, (size_t)n, stdout);
}

void clear_to_eol() {
    fwrite("\033[K", 1, 3, stdout);
}

void flush_input() {
    tcflush(STDIN_FILENO, TCIFLUSH);
}

int read_key() {
    unsigned char c = 0;
    ssize_t n = read(STDIN_FILENO, &c, 1);
    if (n <= 0) return 0;

    if (c != 0x1B) return (int)c;

    fd_set fds;
    FD_ZERO(&fds);
    FD_SET(STDIN_FILENO, &fds);
    struct timeval tv = {0, 50000};
    if (select(STDIN_FILENO + 1, &fds, nullptr, nullptr, &tv) <= 0) {
        return 0x1B;
    }

    unsigned char c2 = 0;
    if (read(STDIN_FILENO, &c2, 1) <= 0) return 0x1B;
    if (c2 != '[') return 0;

    unsigned char c3 = 0;
    if (read(STDIN_FILENO, &c3, 1) <= 0) return 0;

    if (c3 == '<') {
        char numstr[32];
        int  nlen = 0;
        while (nlen < (int)(sizeof(numstr) - 1)) {
            unsigned char ch = 0;
            if (read(STDIN_FILENO, &ch, 1) <= 0) break;
            if (ch == 'M' || ch == 'm') {
                numstr[nlen] = '\0';
                int btn = 0, x = 0, y = 0;
                sscanf(numstr, "%d;%d;%d", &btn, &x, &y);
                g_mouse_x       = x - 1;
                g_mouse_y       = y - 1;
                g_mouse_pressed = (ch == 'M');
                g_mouse_moved   = (btn & 32) != 0;
                return KEY_MOUSE;
            }
            numstr[nlen++] = (char)ch;
        }
        return 0;
    }

    switch (c3) {
        case 'A': return KEY_UP;
        case 'B': return KEY_DOWN;
        case 'C': return KEY_RIGHT;
        case 'D': return KEY_LEFT;
    }

    {
        struct timeval drain_tv = {0, 0};
        fd_set dfds;
        FD_ZERO(&dfds);
        FD_SET(STDIN_FILENO, &dfds);
        while (select(STDIN_FILENO + 1, &dfds, nullptr, nullptr, &drain_tv) > 0) {
            unsigned char dummy;
            if (read(STDIN_FILENO, &dummy, 1) <= 0) break;
            if (dummy >= 0x40 && dummy <= 0x7E) break;
            FD_ZERO(&dfds);
            FD_SET(STDIN_FILENO, &dfds);
            drain_tv = {0, 0};
        }
    }

    return 0;
}

bool read_mouse(MouseEvent& ev) {
    ev.x       = g_mouse_x;
    ev.y       = g_mouse_y;
    ev.pressed = g_mouse_pressed;
    ev.moved   = g_mouse_moved;
    return true;
}

}
