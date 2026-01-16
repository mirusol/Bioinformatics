#A square matrix of an arbitrary size and the corresponding initial vector are given. Implement a software application that makes a prediction on a total of 5 discrete steps, using this matrix and the corresponding vector.


import numpy as np

class StatePredictor:
    def __init__(self, transition_matrix, initial_vector):
        self.matrix = np.array(transition_matrix)
        self.current_state = np.array(initial_vector)
        
        # Validation
        if self.matrix.shape[0] != self.matrix.shape[1]:
            raise ValueError("Input must be a square matrix.")
        if self.matrix.shape[1] != len(self.current_state):
            raise ValueError("Vector size must match matrix dimensions.")

    def predict(self, steps=5):
        predictions = []
        print(f"Initial State (t=0): {self.current_state}")
        
        for t in range(1, steps + 1):
            next_state = np.dot(self.matrix, self.current_state)
            
            predictions.append(next_state)
            self.current_state = next_state
            
            print(f"Step {t}: {self.current_state}")
            
        return predictions


if __name__ == "__main__":
    input_matrix = [
        [0.3, 0.35, 0.35], 
        [0, 0, 1], 
        [0.9, 0, 0.1]
    ]

    initial_v = [1.0, 0.0, 0.0]

    try:
        predictor = StatePredictor(input_matrix, initial_v)
        
        results = predictor.predict(steps=5)
        
    except ValueError as e:
        print(f"Error: {e}")