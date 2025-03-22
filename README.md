# Computer Graphics Project

![GitHub](https://img.shields.io/badge/license-MIT-blue)
![GitHub](https://img.shields.io/badge/language-C%2B%2B-red)

This is a **Computer Graphics Project** built using **C++** and the **Win32 API**. The application allows users to draw various shapes (lines, circles, ellipses) and perform operations like saving, loading, and clearing the drawing canvas. It also includes advanced features like flood fill, polygon clipping, and curve drawing.

---

## Features

- **Drawing Shapes**:
  - Lines (DDA, Midpoint, Parametric Algorithms)
  - Circles (Direct, Polar, Midpoint Algorithms)
  - Ellipses (Direct, Polar Algorithms)
- **Filling Algorithms**:
  - Flood Fill (Recursive and Non-Recursive)
  - Convex and Non-Convex Polygon Filling
- **Clipping**:
  - Line Clipping (Cohen-Sutherland Algorithm)
  - Point Clipping
  - Polygon Clipping
- **Curves**:
  - Hermite Curve
  - Cardinal Spline Curve
- **Color Selection**:
  - Black, Red, Blue, Green
- **File Operations**:
  - Save drawings to a file
  - Load drawings from a file
- **Clear Canvas**:
  - Clear the entire drawing canvas
- **User Interaction**:
  - Left-click to set the starting point
  - Right-click to set the ending point and draw the shape

---

## Getting Started

### Prerequisites

- **Windows OS**: The application is designed to run on Windows.
- **Compiler**: A C++ compiler that supports Win32 API (e.g., Visual Studio, MinGW).

### Installation

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/ahmedm-sallam/Computer-Graphics-Project.git
   cd Computer-Graphics-Project
   ```

2. **Open the Project**:
   - Open the project in your preferred IDE (e.g., Visual Studio).
   - Ensure the project is set up to use the Win32 API.

3. **Build the Project**:
   - Build the project using your IDE or command-line tools.
   - For Visual Studio, press `Ctrl + Shift + B` to build the solution.

4. **Run the Application**:
   - Run the executable file generated in the build directory.

---

## Usage

1. **Drawing Shapes**:
   - Select a drawing mode (e.g., Line, Circle, Ellipse) from the menu.
   - Left-click to set the starting point.
   - Right-click to set the ending point and draw the shape.

2. **Filling Shapes**:
   - Use the flood fill or polygon filling options to fill shapes with color.

3. **Clipping**:
   - Use the clipping options to clip lines, points, or polygons within a specified window.

4. **Changing Colors**:
   - Use the menu to select a color (Black, Red, Blue, Green).

5. **Saving and Loading**:
   - Use the "Save" option to save the current drawing to a file.
   - Use the "Load" option to load a previously saved drawing.

6. **Clearing the Canvas**:
   - Use the "Clear" option to clear the drawing canvas.

---

## Code Structure

- **`main.cpp`**: Contains the main application logic, including the window procedure and drawing functions.
- **Drawing Functions**:
  - `DrawLine_DDA`: Draws a line using the DDA algorithm.
  - `DrawCircle_MidPoint`: Draws a circle using the midpoint algorithm.
  - `DrawEllipse_Direct`: Draws an ellipse using the direct algorithm.
- **Filling Algorithms**:
  - `FloodFill`: Implements the recursive flood fill algorithm.
  - `NRFloodFill`: Implements the non-recursive flood fill algorithm.
- **Clipping Algorithms**:
  - `CohenSuth`: Implements the Cohen-Sutherland line clipping algorithm.
- **File Handling**:
  - `SavePointsToFile`: Saves the current drawing to a file.
  - `LoadPointsFromFile`: Loads a drawing from a file.

---

## Contributing

Contributions are welcome! If you'd like to contribute, please follow these steps:

1. Fork the repository.
2. Create a new branch for your feature or bugfix.
3. Commit your changes.
4. Submit a pull request.

---

## License

This project is licensed under the **MIT License**. See the [LICENSE](LICENSE) file for details.
