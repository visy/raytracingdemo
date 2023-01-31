#pragma once

#include <assert.h>
#include <stdint.h>
#include "core.h"

namespace Pixie
{
    enum MouseButton
    {
        MouseButton_Left	= 0,
        MouseButton_Middle	= 1,
        MouseButton_Right	= 2,
        MouseButton_Num
    };

    enum Key
    {
        Key_Backspace = 0,
        Key_Tab,
        Key_Enter,
        Key_Escape,

        Key_Up,
        Key_Down,
        Key_Left,
        Key_Right,

        Key_Home,
        Key_End,
        Key_PageUp,
        Key_PageDown,

        Key_Delete,
        Key_Insert,

        Key_LeftShift,
        Key_RightShift,
        Key_LeftControl,
        Key_RightControl,
        Key_LeftAlt,
        Key_RightAlt,

        Key_F1,
        Key_F2,
        Key_F3,
        Key_F4,
        Key_F5,
        Key_F6,
        Key_F7,
        Key_F8,
        Key_F9,
        Key_F10,
        Key_F11,
        Key_F12,

        Key_ASCII_Start = 32,

        // 32 to 127 are ASCII printable characters.
        // Note: these are unshifted keys.

        Key_ASCII_End = 127,

        Key_Num
    };

    enum
    {
        MaxPlatformKeys = 256
    };

    class Window
    {
        public:
            Window();
            ~Window();

            // Open the Pixie window with the specified title bar, width, and height.
            // If scale is greater than 1 the window will be rendered scale times larger
            // and the buffer will be stretched to fit.
            bool Open(const char* title, int width, int height, bool fullscreen, bool maintainAspectRatio = false, int scale = 1);

            // Close the Pixie window.
            void Close();

            // Update the Pixie window. This will copy the backing buffer to the actual window.
            bool Update();

            // Returns true in the frame the mouse button went down.
            bool HasMouseGoneDown(MouseButton button) const;

            // Returns true in the frame the mouse button went up.
            bool HasMouseGoneUp(MouseButton button) const;

            // Returns true if the mouse button is currently down.
            bool IsMouseDown(MouseButton button) const;

            // Returns true in the frame the key specified went down.
            bool HasKeyGoneDown(Key key) const;

            // Returns true if any key went down in the frame.
            bool HasAnyKeyGoneDown() const;

            // Returns true in the frame the key specified went up.
            bool HasKeyGoneUp(Key key) const;

            // Returns true if the key specified is currently down.
            bool IsKeyDown(Key key) const;

            // Returns true if any key is currently down.
            bool IsAnyKeyDown() const;

            // Returns an array of ASCII input for the current frame.
            const char* GetInputCharacters() const;

            // Clears the ASCII input for the current frame.
            void ClearInputCharacters();

            // Returns the current mouse X position.
            int GetMouseX() const;

            // Returns the current mouse Y position.
            int GetMouseY() const;

            // Returns the time delta in seconds since the last time the window was updated.
            float GetDelta() const;

            // Returns the time in seconds since the window was opened.
            float GetTime() const;

            // Returns the backing buffer for the window.
            uint32_t* GetPixels() const;

            // Returns the width of the window.
            uint32_t GetWidth() const;

            // Returns the height of the window.
            uint32_t GetHeight() const;

            // Returns the scale of the window.
            uint32_t GetScale() const;

            // Key callback handler. Called on any key state change.
            typedef void(*KeyCallback)(Key key, bool down);
            void SetKeyCallback(KeyCallback callback);

            // Used by the window procedure to update key and mouse state.
            void SetMouseButtonDown(MouseButton button, bool down);
            void SetKeyDown(int key, bool down);
            void AddInputCharacter(char c);

        private:
            void PlatformInit();
            bool PlatformOpen(const char* title, int width, int height);
            bool PlatformUpdate();
            void PlatformClose();

            void UpdateMouse();
            void UpdateKeyboard();

            int m_mouseX;
            int m_mouseY;
            bool m_lastMouseButtonDown[MouseButton_Num];
            bool m_mouseButtonDown[MouseButton_Num];

            int m_keyMap[Key_Num];
            bool m_lastKeyDown[MaxPlatformKeys];
            bool m_keyDown[MaxPlatformKeys];
            bool m_anyKeyDown;
            char m_inputCharacters[16+1];

            float m_delta;

            uint32_t* m_pixels;
            uint32_t m_width;
            uint32_t m_height;
            uint32_t m_windowWidth;
            uint32_t m_windowHeight;
            int m_scale;
            bool m_fullscreen;
            bool m_maintainAspectRatio;
            float m_scalex, m_scaley;

            void* m_window;

            float m_time;
            int64_t m_lastTime;
            int64_t m_freq;

            KeyCallback m_keyCallback;
    };

    inline int Window::GetMouseX() const
    {
        return m_mouseX;
    }

    inline int Window::GetMouseY() const
    {
        return m_mouseY;
    }

    inline float Window::GetDelta() const
    {
        return m_delta;
    }

    inline float Window::GetTime() const
    {
        return m_time;
    }

    inline uint32_t* Window::GetPixels() const
    {
        return m_pixels;
    }

    inline uint32_t Window::GetWidth() const
    {
        return m_width;
    }

    inline uint32_t Window::GetHeight() const
    {
        return m_height;
    }

    inline uint32_t Window::GetScale() const
    {
        return m_scale;
    }

    inline bool Window::HasMouseGoneDown(MouseButton button) const
    {
        return !m_lastMouseButtonDown[button] && m_mouseButtonDown[button];
    }

    inline bool Window::HasMouseGoneUp(MouseButton button) const
    {
        return m_lastMouseButtonDown[button] && !m_mouseButtonDown[button];
    }

    inline bool Window::IsMouseDown(MouseButton button) const
    {
        return m_mouseButtonDown[button];
    }

    inline bool Window::HasAnyKeyGoneDown() const
    {
        for (int i = 0; i < Key_Num; i++)
        {
            if (HasKeyGoneDown((Key)i))
                return true;
        }

        return false;
    }

    inline bool Window::HasKeyGoneDown(Key key) const
    {
        int index = m_keyMap[key];
        if (index == -1)
            return false;
        assert(index >= 0 && index < MaxPlatformKeys);
        return !m_lastKeyDown[index] && m_keyDown[index];
    }

    inline bool Window::HasKeyGoneUp(Key key) const
    {
        int index = m_keyMap[key];
        if (index == -1)
            return false;
        assert(index >= 0 && index < MaxPlatformKeys);
        return m_lastKeyDown[index] && !m_keyDown[index];
    }

    inline bool Window::IsKeyDown(Key key) const
    {
        int index = m_keyMap[key];
        if (index == -1)
            return false;
        assert(index >= 0 && index < MaxPlatformKeys);
        return m_keyDown[index];
    }

    inline bool Window::IsAnyKeyDown() const
    {
        for (int i = 0; i < Key_Num; i++)
        {
            if (IsKeyDown((Key)i))
                return true;
        }

        return false;
    }

    inline const char* Window::GetInputCharacters() const
    {
        return m_inputCharacters;
    }

    inline void Window::ClearInputCharacters()
    {
        m_inputCharacters[0] = 0;
    }

    inline void Window::SetMouseButtonDown(MouseButton button, bool down)
    {
        m_mouseButtonDown[button] = down;
    }

    inline void Window::SetKeyDown(int platformKey, bool down)
    {
        assert(platformKey >= 0 && platformKey < MaxPlatformKeys);
        if (m_keyDown[platformKey] == down)
            return;

        m_keyDown[platformKey] = down;
        if (m_keyCallback)
        {
            for (int i = 0; i < Key_Num; i++)
            {
                int key = m_keyMap[i];
                if (key == platformKey)
                {
                    m_keyCallback((Key)i, down);
                    return;
                }
            }
        }
    }

    inline void Window::SetKeyCallback(KeyCallback callback)
    {
        m_keyCallback = callback;
    }
}
