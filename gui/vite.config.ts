import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'

// Tauri drives the dev server on a fixed port and needs a stable host.
export default defineConfig({
  plugins: [react()],
  clearScreen: false,
  server: {
    port: 1420,
    strictPort: true,
    // The license text is imported from the repository root (../LICENSE).
    fs: { allow: ['..'] },
  },
  build: {
    target: 'safari14',
    sourcemap: true,
  },
})
