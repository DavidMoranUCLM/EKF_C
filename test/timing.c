#include "EKF.h"
#include "esp_log.h"
#include "tictoc.h"
#include "unity.h"

EKF_ctx_t ekf = {0};
measures_t meas = {0};
tictoc_t *initTic;
tictoc_t *stepTic;

const char *TAG = "EKF_C";

TEST_CASE("1", "[EKF_C]") {
  ESP_LOGI(TAG,"Beginning test...");
  for (int i = 0; i < 100; i++) {
    //ESP_LOGI(TAG,"Iteration %d", i);
    tic(initTic);
    ekfInit(&ekf, &meas);
    toc(initTic);
    ekfDeinit(&ekf);
    ESP_LOGI(TAG, "Free heap size: %lu B", esp_get_free_heap_size());
  }
  tictoc_print(initTic);
}

void setUp(void) {
  meas.acc[0] = 0;
  meas.acc[1] = 0;
  meas.acc[2] = 9.81;
  meas.mag[0] = 0;
  meas.mag[1] = 1;
  meas.mag[2] = 0;
  initTic = tictoc_new("Init time");
  stepTic = tictoc_new("Step time");
  ESP_LOGI(TAG,"Timers Set");
}