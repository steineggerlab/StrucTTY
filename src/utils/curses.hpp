#pragma once

#if __has_include(<ncurses.h>)
  #include <ncurses.h>
#elif __has_include(<curses.h>)
  #include <curses.h>
#elif __has_include(<ncurses/ncurses.h>)
  #include <ncurses/ncurses.h>
#elif __has_include(<ncurses/curses.h>)
  #include <ncurses/curses.h>
#elif __has_include(<ncursesw/ncurses.h>)
  #include <ncursesw/ncurses.h>
#elif __has_include(<ncursesw/curses.h>)
  #include <ncursesw/curses.h>
#else
  #error "No curses/ncurses header found. Install ncurses development package."
#endif
